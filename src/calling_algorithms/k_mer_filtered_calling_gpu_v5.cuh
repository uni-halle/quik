//
// Created by steffen on 03.02.26.
//

#ifndef INC_2OPT_KMER_FILTERED_CALLING_GPU_V5_CUH
#define INC_2OPT_KMER_FILTERED_CALLING_GPU_V5_CUH

#include "barcode_calling_algorithm.h"
#include "../distance/weighted_levenshtein_v1.cuh"
#include "../distance/weighted_sequence_levenshtein_v1.cuh"
#include "../cuda_helper.cuh"

namespace barcode_calling {

    namespace k_mer_filtered_calling_gpu_v5_detail {

        template <unsigned CANDIDATE_COUNT>
        inline size_t compute_shared_bytes(unsigned barcodes_per_chunk) {

            auto align_up = [](size_t x, size_t a) {
                return (x + a - 1) & ~(a - 1);
            };

            size_t offset = 0;

            // pseudo_distance
            offset = align_up(offset, alignof(int32_t));
            offset += barcodes_per_chunk * sizeof(int32_t);

            // candidate_barcodes
            offset = align_up(offset, alignof(unsigned));
            offset += (CANDIDATE_COUNT + 1) * sizeof(unsigned);

            // candidate_distance
            offset = align_up(offset, alignof(int32_t));
            offset += (CANDIDATE_COUNT + 1) * sizeof(int32_t);

            // barcode_ids
            offset = align_up(offset, alignof(unsigned));
            offset += CANDIDATE_COUNT * sizeof(unsigned);

            return offset;
        }


        /**
         * Determine from the set of indices an index i such that keys[i] is minimum.
         *
         * @param index_set Subset of indices. (Will be overridden)
         * @param keys Set of keys.
         * @param size Size of index set.
         * @return
         */
        __device__ unsigned arg_min_intra_block(unsigned* index_set, const int32_t* keys, unsigned size) {
            unsigned active_elements = size;
            while (active_elements > 1) {
                unsigned half = (active_elements + 1) / 2; // half = ceil(active_elements/2)
                for (unsigned i = threadIdx.x; half + i < active_elements; i += warpSize) {
                    unsigned j = half + i;
                    if (keys[index_set[i]] > keys[index_set[j]])
                        index_set[i] = index_set[j];
                }
                active_elements = half;
            }
            return index_set[0];
        }

        template <typename distance_function,
                  unsigned k, // length of k-mers
                  unsigned barcodes_per_chunk,
                  unsigned PSEUDO_DISTANCE_WINDOW_SIZE = DEFAULT_PSEUDO_DISTANCE_WINDOW_SIZE,
                  unsigned CANDIDATE_COUNT = DEFAULT_CANDIDATE_COUNT
        >
        __global__ void kernel(

            // input variables
            const barcode* barcodes,
            const size_t* read_start_index,
            const uint8_t* read_data,
            const unsigned* read_lengths,
            unsigned** k_mer_map_entries, // k-mer data structures (chunk_count many)
            size_t** k_mer_map_start_indices, // k-mer data structures (chunk_count many)
            alignment_costs costs,
            const size_t barcode_count, // number of barcodes
            const size_t read_count, // number of reads

            // output variables
            unsigned* out_barcodes, // out: 2*read_count entries
            int32_t* out_distances // out: 2*read_count entries
        ) {
            const unsigned chunk_count = (barcode_count + barcodes_per_chunk - 1) / barcodes_per_chunk;

            const uint32_t k_mer_count = k_mer_map_host_v1<k>::get_k_mer_count(); // 4^k
            const uint32_t mask = (1 << 2 * k) - 1; // bitmask with ones at its 2k least significant positions

            assert(CANDIDATE_COUNT <= barcode_count);
            assert(barcodes_per_chunk >= CANDIDATE_COUNT);

            /******************************************************************************************
             * Each block is associated to a single read r.
             *
             * In the first step, the block determines a list of barcodes with smallest pseudo distance.
             * These barcodes are candidates for the final assignment.
             *****************************************************************************************/

            const unsigned read_id = blockIdx.x;

            // construct read from global memory
            const uint8_t* read_data_ptr = read_data + read_start_index[read_id];
            const read r(read_data_ptr, read_lengths[read_id]);

            /**********************************************************************************
             * Prepare the shared memory.
             *
             * 1. The first part is reserved for the array of pseudo distances.
             *
             *        pseudo_distance = int32_t[barcodes_per_chunk]
             *
             * 2. The second part contains the candidate barcodes with associated distances.
             *
             *       candidate_barcodes = unsigned[CANDIDATE_COUNT+1]
             *       candiate_distances = int32_t[CANDIDATE_COUNT+1]
             *
             * 3. The third part is used as temporary memory for reduction operations.
             *
             *       barcode_ids = unsigned[CANDIDATE_COUNT]
             *********************************************************************************/

            extern __shared__ unsigned char smem[];

            // alignment helper
            auto align_up = [](size_t x, size_t a) {
                return (x + a - 1) & ~(a - 1);
            };

            size_t offset = 0;

            // pseudo_distance
            offset = align_up(offset, alignof(int32_t));
            int32_t* pseudo_distance = reinterpret_cast<int32_t*>(smem + offset);
            offset += barcodes_per_chunk * sizeof(int32_t);

            // candidate_barcodes
            offset = align_up(offset, alignof(unsigned));
            unsigned* candidate_barcodes = reinterpret_cast<unsigned*>(smem + offset);
            offset += (CANDIDATE_COUNT + 1) * sizeof(unsigned);

            // candidate_distance
            offset = align_up(offset, alignof(int32_t));
            int32_t* candidate_distance = reinterpret_cast<int32_t*>(smem + offset);
            offset += (CANDIDATE_COUNT + 1) * sizeof(int32_t);

            // barcode_ids
            offset = align_up(offset, alignof(unsigned));
            unsigned* barcode_ids = reinterpret_cast<unsigned*>(smem + offset);
            offset += CANDIDATE_COUNT * sizeof(unsigned);

            /****************************************************************************
             * initialize the candidate list with dummy elements
             ***************************************************************************/

            for (unsigned i = threadIdx.x; i <= CANDIDATE_COUNT; i += warpSize) {
                candidate_barcodes[i] = UINT_MAX;
                candidate_distance[i] = INT32_MAX;
            }

            /****************************************************************************
             * For each chunk of barcodes
             ***************************************************************************/

            for (unsigned chunk_id = 0; chunk_id < chunk_count; ++chunk_id) {

                unsigned long long t2 = clock64();

                /****************************************************************************
                 * Step 1: Compute the pseudo distance to all barcodes of the current chunk.
                 ***************************************************************************/

                // initialize the pseudo distance array
                for (unsigned barcode_id = threadIdx.x; barcode_id < barcodes_per_chunk; barcode_id += warpSize) {
                    pseudo_distance[barcode_id] = 0;
                }

                // prepare the first k-1 characters of k-mer m = r[0...k-1]
                uint32_t m = 0; // m = r[0...k-2]
                for (unsigned l = 0; l < k - 1; l++)
                    m = (m << 2) + r[l].to_uint8();

                // for each position i in which a k-mer may start in the read
                for (int i = 0; i + k <= r.length(); i++) {

                    // invariant: m == r[i...i+k-1]
                    m = ((m << 2) & mask) + r[i + k - 1].to_uint8();

                    // for each position j in a window around position i
                    int j = 0;
                    if (i > PSEUDO_DISTANCE_WINDOW_SIZE)
                        j = i - PSEUDO_DISTANCE_WINDOW_SIZE;
                    for (; j + k <= BARCODE_LENGTH && j <= i + PSEUDO_DISTANCE_WINDOW_SIZE; j++) {

                        size_t start = k_mer_map_start_indices[chunk_id][j * k_mer_count + m];
                        size_t stop = k_mer_map_start_indices[chunk_id][j * k_mer_count + m + 1];

                        // for each barcodes in which m occurs at position j
                        for (size_t l = start + threadIdx.x; l < stop; l += warpSize) {
                            unsigned barcode_id = k_mer_map_entries[chunk_id][l];
                            assert(barcode_id < barcodes_per_chunk);
                            // no atomic operation needed as barcodes do not occur twice in the same list
                            pseudo_distance[barcode_id] += abs(i - j) - BARCODE_LENGTH;
                        }
                    }
                }

                /*****************************************************************************************
                 * Step 2: Use the newly computed pseudo distances into update the candidate list.
                 *
                 * For this purpose, we keep the candidate list sorted by their pseudo distance.
                 * Iteratively, we determine minimum barcodes from the current chunk and insert them
                 * into the candidate list.
                 ****************************************************************************************/

                for (unsigned rank = 0; rank < CANDIDATE_COUNT; rank++) {

                    /*****************************************************************************************
                     * Partition the barcodes of this chunk between the threads.
                     * Each thread determines a barcode with minimum pseudo within its partition.
                     *
                     * This is mostly redundant work as it has to repeated for each rank.
                     * However, it seems to be faster than each thread storing a sorted list of barcode ids.
                     * (Which is what v6 does.)
                     ****************************************************************************************/

                    unsigned my_best_barcode_id = threadIdx.x;
                    for (unsigned barcode_id = threadIdx.x + warpSize;
                         barcode_id < barcodes_per_chunk; barcode_id += warpSize) {
                        if (pseudo_distance[my_best_barcode_id] > pseudo_distance[barcode_id])
                            my_best_barcode_id = barcode_id;
                    }

                    /*****************************************************************************************
                     * Find a barcode with minimum pseudo distance within the current chunk.
                     ****************************************************************************************/

                    barcode_ids[threadIdx.x] = my_best_barcode_id;
                    unsigned best_barcode_in_chunk = arg_min_intra_block(barcode_ids, pseudo_distance, warpSize);
                    assert(best_barcode_in_chunk < barcodes_per_chunk);
                    int32_t best_pseudo_dist_in_chunk = pseudo_distance[best_barcode_in_chunk];

                    /*****************************************************************************************
                    * Test whether our newly found barcode is good enough to fit into the candidate list.
                    *
                    * For this purpose, we need to test whether the worst candidate barcode has a larger
                    * pseudo distance then the best in the current chunk.
                    ****************************************************************************************/

                    if (best_pseudo_dist_in_chunk >= candidate_distance[CANDIDATE_COUNT - 1])
                        break;

                    /*****************************************************************************************
                     * Set this pseudo distance to infinity so it cannot be used a second time.
                     ****************************************************************************************/

                    if (threadIdx.x == 0)
                        pseudo_distance[best_barcode_in_chunk] = INT32_MAX;
                    __threadfence_block();

                    /*****************************************************************************************
                     * Thread 0 insertion sorts the new barcode into the shared list of candidate barcodes.
                     ****************************************************************************************/

                    if (threadIdx.x == 0) {

                        unsigned j = CANDIDATE_COUNT;
                        for (; j > 0 && candidate_distance[j - 1] > best_pseudo_dist_in_chunk; --j) {
                            candidate_barcodes[j] = candidate_barcodes[j - 1];
                            candidate_distance[j] = candidate_distance[j - 1];
                        }

                        unsigned global_barcode_id = best_barcode_in_chunk + chunk_id * barcodes_per_chunk;
                        candidate_barcodes[j] = global_barcode_id;
                        candidate_distance[j] = best_pseudo_dist_in_chunk;
                    }
                }
            }

            /****************************************************************************************************
             * Step 3: Replace the pseudo distances by the exact distance.
             ***************************************************************************************************/

            for (unsigned i = threadIdx.x; i < CANDIDATE_COUNT; i += warpSize) {
                barcode b = barcodes[candidate_barcodes[i]];
                candidate_distance[i] = distance_function::evaluate(b, r, costs);
            }

            /*********************************************************************************
             * Step 4: Determine the barcode with minimum exact distance.
             *********************************************************************************/

            for (unsigned i = threadIdx.x; i < CANDIDATE_COUNT; i += warpSize)
                barcode_ids[i] = i;

            // in-block reduction
            unsigned i = arg_min_intra_block(barcode_ids, candidate_distance, CANDIDATE_COUNT);

            if (threadIdx.x == 0) {

                unsigned barcode_id_1 = candidate_barcodes[i];
                int32_t distance_1 = candidate_distance[i];

                // write to global memory
                out_barcodes[read_id] = barcode_id_1;
                out_distances[read_id] = distance_1;

                // make barcode_id_1 invalid
                candidate_distance[i] = INT32_MAX;
            }

            /*********************************************************************************
             * Step 5: Determine the barcode with the second to minimum exact distance.
             *********************************************************************************/

            for (unsigned i = threadIdx.x; i < CANDIDATE_COUNT; i += warpSize)
                barcode_ids[i] = i;

            // in-block reduction
            unsigned j = arg_min_intra_block(barcode_ids, candidate_distance, CANDIDATE_COUNT);

            if (threadIdx.x == 0) {

                unsigned barcode_id_2 = candidate_barcodes[j];
                int32_t distance_2 = candidate_distance[j];

                // write to global memory
                out_barcodes[read_id + read_count] = barcode_id_2;
                out_distances[read_id + read_count] = distance_2;
            }
        }
    }


    template <unsigned k,
              unsigned barcodes_per_chunk = 2048>
    class k_mer_filtered_calling_gpu_v5 : public barcode_calling_algorithm {

    public:
        explicit k_mer_filtered_calling_gpu_v5(const distance_measure& dist, unsigned rejection_threshold)
            : barcode_calling_algorithm(std::to_string(k) + "_mer_filtered_calling_gpu_v5<"
                                        + std::to_string(barcodes_per_chunk) + ">",
                                        dist, rejection_threshold) {
            if (!dist.get_name().starts_with("weighted"))
                throw std::runtime_error("Unsupported distance measure!");
        }

        /**
         * Run the algorithm
         **/
        extended_barcode_assignment run(
            const barcode_set& barcodes,
            const read_set& reads
        ) const override {

            if (dist.get_name() == "weighted_levenshtein_v1")
                return run<weighted_levenshtein_v1>(barcodes, reads, unit_costs());
            if (dist.get_name() == "weighted_sequence_levenshtein_v1")
                return run<weighted_sequence_levenshtein_v1>(barcodes, reads, unit_costs());

            throw std::runtime_error("Unsupported distance measure!");
        }

        /**
         *  Do the calculation.
         **/
        template <typename distance_measure,
                  unsigned PSEUDO_DISTANCE_WINDOW_SIZE = DEFAULT_PSEUDO_DISTANCE_WINDOW_SIZE,
                  unsigned CANDIDATE_COUNT = DEFAULT_CANDIDATE_COUNT
        >
        extended_barcode_assignment run(
            const barcode_set& barcodes,
            const read_set& reads,
            const alignment_costs& costs) const {

            extended_barcode_assignment ass(reads.size()); // will be returned

            /***********************************************************************
             * The main idea of this algorithm is to store the pseudo distances
             * within the GPU's shared memory instead of the global memory.
             * In doing so, the calculation of the pseudo distances may be faster.
             *
             * For this purpose, we divide the barcode set into chunks of
             * approximately the same size. The chunks are processed sequentially.
             *
             * The shared memory of each block is organized as follows:
             *
             *   1. The first part is reserved for the array of pseudo distances.
             *
             *        pseudo_distance = int32_t[barcodes_per_chunk]
             *
             * 2. The second part contains the candidate barcodes with associated distances.
             *
             *       candidate_barcodes = unsigned[CANDIDATE_COUNT+1]
             *       candiate_distances = int32_t[CANDIDATE_COUNT+1]
             *
             * 3. The third part is used as temporary memory for reduction operations.
             *
             *       barcode_ids = unsigned[CANDIDATE_COUNT]
             ************************************************************************/

            size_t shared_size = k_mer_filtered_calling_gpu_v5_detail::
                compute_shared_bytes<CANDIDATE_COUNT>(barcodes_per_chunk);

            // determine the size of the shared memory and the number of barcodes chunks
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, 0);

            /*std::cout << "required shared memory per block:  " << shared_size << std::endl;
            std::cout << "available shared memory per block: " << prop.sharedMemPerBlock << std::endl;*/

            assert(shared_size <= prop.sharedMemPerBlock);

            /*************************************************************************
             * Partition the barcode set into chunks
             ************************************************************************/

            std::vector<barcode_set> barcode_chunks;
            for (unsigned barcode_start_id = 0; barcode_start_id < barcodes.size();
                 barcode_start_id += barcodes_per_chunk) {

                /*************************************************************************
                 * Create a new chunk consisting of the original barcodes
                 * barcodes[barcode_start_id ... barcode_start_id + barcodes_per_chunk-1]
                 ************************************************************************/

                barcode_set new_chunk;
                for (unsigned barcode_id = barcode_start_id;
                     barcode_id < barcode_start_id + barcodes_per_chunk && barcode_id < barcodes.size();
                     barcode_id++) {
                    new_chunk.add(barcodes[barcode_id], barcodes.get_name_of(barcode_id));
                }
                barcode_chunks.push_back(new_chunk);
            }

            /***********************************************************************
             * Allocate gpu memory for the barcodes, reads, etc.
             ************************************************************************/

            // gpu input data
            barcode* barcodes_dev;
            size_t* read_start_index_dev;
            uint8_t* read_sequence_data_dev;
            unsigned* read_lengths_dev;

            CUDA_CHECK(cudaMalloc(&barcodes_dev, barcodes.size() * sizeof(barcode)));
            CUDA_CHECK(cudaMalloc(&read_start_index_dev, reads.get_start_indices().size() * sizeof(size_t)));
            CUDA_CHECK(cudaMalloc(&read_sequence_data_dev, reads.get_sequence_data().size() * sizeof(uint8_t)));
            CUDA_CHECK(cudaMalloc(&read_lengths_dev, reads.get_sequence_lengths().size() * sizeof(unsigned)));

            /************************************************************************************************
            * Copy the barcodes and reads to gpu memory.
            ************************************************************************************************/

            // barcodes
            CUDA_CHECK(cudaMemcpy(barcodes_dev, barcodes.data(),
                barcodes.size() * sizeof(barcode), cudaMemcpyHostToDevice));

            // read start positions
            CUDA_CHECK(cudaMemcpy(read_start_index_dev, reads.get_start_indices().data(),
                reads.get_start_indices().size() * sizeof(size_t),
                cudaMemcpyHostToDevice));

            // concatenated read lengths
            CUDA_CHECK(cudaMemcpy(read_lengths_dev, reads.get_sequence_lengths().data(),
                reads.get_sequence_lengths().size() * sizeof(unsigned),
                cudaMemcpyHostToDevice));

            // concatenated read sequence data
            CUDA_CHECK(cudaMemcpy(read_sequence_data_dev, reads.get_sequence_data().data(),
                reads.get_sequence_data().size() * sizeof(uint8_t),
                cudaMemcpyHostToDevice));

            /************************************************************************************************
             * Allocate the output memory for the calling kernel.
             *
             * For each read r, the kernel outputs the 2 closest barcodes together with the distances.
             ************************************************************************************************/

            // gpu output data
            unsigned* out_barcodes_dev;
            int32_t* out_distances_dev;

            CUDA_CHECK(cudaMalloc(&out_barcodes_dev, 2*reads.size() * sizeof(unsigned)));
            CUDA_CHECK(cudaMalloc(&out_distances_dev, 2*reads.size() * sizeof(int32_t)));

            /****************************************************************************
             * Prepare the k-mer data structures of each chunk seperately.
             *
             * For each combination of a k-mer x in {A,C,G,T}^k and position i
             * in [0...BARCODE_LENGTH-1], we construct a list of barcodes b in which
             * x occurs in b at position i.
             *
             * Thus, k_mer_map[x][i] is a list of barcodes b for which x == b[i...i+k-1].
             *
             * We need to copy the k-mer maps to GPU memory.
             *
             * For this purpose, we need to change the data layout.
             ****************************************************************************/

            // for each barcode chunk a pointer to the associated k-mer map in GPU memory
            std::vector<unsigned*> k_mer_map_entries_host;
            std::vector<size_t*> k_mer_map_start_indices_host;

            for (const barcode_set& barcode_chunk : barcode_chunks) {

                /***********************************************************************
                 * Create the k_mer map for the current chunk in host memory.
                 *
                 * For the j-th k-mer map (associated to the j-th chunk), we concatenate the
                 * individual barcode lists k_mer_map[j][x][i] into a single array
                 * k_mer_map_entries_dev[j], similarly to an adjacency array of a graph.
                 *
                 * The barcode indices from k_map_map[j][x][i] are then located in the
                 * array k_mer_map_entries_dev[j] at the positions
                 *
                 *  k_mer_map_start_index_dev[j][i * k_mer_count + x] to
                 *  k_mer_map_start_index_dev[j][i * k_mer_count + x + 1] - 1.
                 ****************************************************************************/

                k_mer_map_host_v1<k> k_mer_map(barcode_chunk);
                const size_t k_mer_count = k_mer_map.get_k_mer_count();

                const unsigned k_mer_start_position_count = BARCODE_LENGTH - k + 1;
                std::vector<size_t> current_k_mer_map_start_index_host(k_mer_count * k_mer_start_position_count + 1);
                std::vector<unsigned> current_k_mer_map_entries_host(barcode_chunk.size() * k_mer_start_position_count);

                size_t running_index = 0;

                // for each k-mer start position i within the barcode
                for (unsigned i = 0; i <= BARCODE_LENGTH - k; i++) {

                    // for each possible k-mer x
                    for (unsigned x = 0; x < k_mer_count; x++) {
                        // concat the lists into a single array
                        memcpy(&current_k_mer_map_entries_host[running_index],
                               k_mer_map[x][i].data(),
                               k_mer_map[x][i].size() * sizeof(unsigned));
                        running_index += k_mer_map[x][i].size();
                        current_k_mer_map_start_index_host[i * k_mer_count + x + 1] = running_index;
                    }
                }

                /***************************************************************************
                 * Transfer the k-mer map to the GPU memory
                 **************************************************************************/

                unsigned* current_k_mer_map_entries_dev;
                size_t* current_k_mer_map_start_index_dev;

                // allocate GPU memory
                CUDA_CHECK(cudaMalloc(&current_k_mer_map_entries_dev,
                    current_k_mer_map_entries_host.size() * sizeof(unsigned)));
                CUDA_CHECK(cudaMalloc(&current_k_mer_map_start_index_dev,
                    current_k_mer_map_start_index_host.size() * sizeof(size_t)));

                // copy the data to the GPU memory
                CUDA_CHECK(cudaMemcpy(current_k_mer_map_entries_dev, current_k_mer_map_entries_host.data(),
                    current_k_mer_map_entries_host.size() * sizeof(unsigned),
                    cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemcpy(current_k_mer_map_start_index_dev, current_k_mer_map_start_index_host.data(),
                    current_k_mer_map_start_index_host.size() * sizeof(size_t),
                    cudaMemcpyHostToDevice));

                // store the pointers to the GPU memory
                k_mer_map_entries_host.push_back(current_k_mer_map_entries_dev);
                k_mer_map_start_indices_host.push_back(current_k_mer_map_start_index_dev);
            }

            /***************************************************************************
             * Copy the pointers to the k-mer maps to the GPU
             **************************************************************************/

            unsigned** k_mer_map_entries_dev;
            size_t** k_mer_map_start_indices_dev;

            CUDA_CHECK(cudaMalloc(&k_mer_map_entries_dev,
                k_mer_map_entries_host.size() * sizeof(unsigned*)));
            CUDA_CHECK(cudaMalloc(&k_mer_map_start_indices_dev,
                k_mer_map_start_indices_host.size() * sizeof(size_t*)));

            CUDA_CHECK(cudaMemcpy(k_mer_map_entries_dev, k_mer_map_entries_host.data(),
                k_mer_map_entries_host.size() * sizeof(unsigned*), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(k_mer_map_start_indices_dev, k_mer_map_start_indices_host.data(),
                k_mer_map_start_indices_host.size() * sizeof(size_t*), cudaMemcpyHostToDevice));

            /************************************************************************************************
             * Call the calling kernel.
             ************************************************************************************************/

            k_mer_filtered_calling_gpu_v5_detail::kernel
                <distance_measure, k, barcodes_per_chunk, PSEUDO_DISTANCE_WINDOW_SIZE, CANDIDATE_COUNT>
                <<<reads.size(), 32, shared_size>>>
                (barcodes_dev, read_start_index_dev, read_sequence_data_dev, read_lengths_dev,
                 k_mer_map_entries_dev, k_mer_map_start_indices_dev,
                 costs, barcodes.size(), reads.size(),
                 out_barcodes_dev, out_distances_dev);

            CUDA_CHECK(cudaDeviceSynchronize());

            /************************************************************************************************
             * Copy the results to the host memory.
             ************************************************************************************************/

            std::vector<unsigned> barcodes_host(reads.size());
            std::vector<int32_t> distances_host(reads.size());

            // closest barcodes with distances
            CUDA_CHECK(cudaMemcpy(barcodes_host.data(), out_barcodes_dev,
                reads.size() * sizeof(unsigned), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(distances_host.data(), out_distances_dev,
                reads.size() * sizeof(int32_t), cudaMemcpyDeviceToHost));

            for (unsigned read_id = 0; read_id < reads.size(); read_id++) {
                if (distances_host[read_id] <= rejection_threshold)
                    ass.assign_as_1st_barcode(read_id, barcodes_host[read_id], distances_host[read_id]);
            }

            // 2nd closest barcodes with distances
            CUDA_CHECK(cudaMemcpy(barcodes_host.data(), out_barcodes_dev + reads.size(),
                reads.size() * sizeof(unsigned), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(distances_host.data(), out_distances_dev + reads.size(),
                reads.size() * sizeof(int32_t), cudaMemcpyDeviceToHost));

            for (unsigned read_id = 0; read_id < reads.size(); read_id++) {
                if (distances_host[read_id] <= rejection_threshold)
                    ass.assign_as_2nd_barcode(read_id, barcodes_host[read_id], distances_host[read_id]);
            }

            /************************************************************************************************
             * Free the memory.
             ************************************************************************************************/

            CUDA_CHECK(cudaFree(barcodes_dev));
            CUDA_CHECK(cudaFree(read_start_index_dev));
            CUDA_CHECK(cudaFree(read_sequence_data_dev));
            CUDA_CHECK(cudaFree(read_lengths_dev));
            CUDA_CHECK(cudaFree(out_barcodes_dev));
            CUDA_CHECK(cudaFree(out_distances_dev));
            CUDA_CHECK(cudaFree(k_mer_map_entries_dev));
            CUDA_CHECK(cudaFree(k_mer_map_start_indices_dev));
            for (unsigned chunk_id = 0; chunk_id < barcode_chunks.size(); chunk_id++) {
                CUDA_CHECK(cudaFree(k_mer_map_entries_host[chunk_id]));
                CUDA_CHECK(cudaFree(k_mer_map_start_indices_host[chunk_id]));
            }

            return ass;
        }


    };


}


#endif //INC_2OPT_KMER_FILTERED_CALLING_GPU_V5_CUH
