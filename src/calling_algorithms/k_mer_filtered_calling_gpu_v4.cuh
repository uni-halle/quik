//
// Created by steffen on 25.06.24.
//

#ifndef INC_2OPT_KMER_FILTERED_CALLING_GPU_V4_CUH
#define INC_2OPT_KMER_FILTERED_CALLING_GPU_V4_CUH

#include "barcode_calling_algorithm.h"
#include "../distance/weighted_levenshtein_v1.cuh"
#include "../distance/weighted_sequence_levenshtein_v1.cuh"
#include "../cuda_helper.cuh"

namespace barcode_calling {

    namespace k_mer_filtered_calling_gpu_v4_detail {

        template <typename distance_function,
                  unsigned k, // length of k-mers
                  unsigned threads_per_block = 256,
                  unsigned PSEUDO_DISTANCE_WINDOW_SIZE = DEFAULT_PSEUDO_DISTANCE_WINDOW_SIZE,
                  unsigned MAX_INDEX_SIZE = DEFAULT_CANDIDATE_COUNT
        >
        __global__ void kernel(

            // input variables
            const barcode* barcodes,
            const size_t* read_start_index,
            const uint8_t* read_data,
            const unsigned* read_lengths,
            const unsigned* k_mer_map_entries, // k-mer data structure
            const unsigned* k_mer_start_index, // k-mer data structure
            int32_t* pseudo_distances, // for each thread: array of size barcode_count
            unsigned* barcode_candidates, // for each thread: array of size barcode_count
            alignment_costs costs,
            const size_t barcode_count, // number of barcodes
            const size_t read_count, // number of reads

            // output variables
            unsigned* out_barcodes, // out: 2*read_count entries
            int32_t* out_distances // out: 2*read_count entries
        ) {
            const uint32_t k_mer_count = k_mer_map_host_v1<k>::get_k_mer_count(); // 4^k
            const uint32_t mask = (1 << 2 * k) - 1;
            // bitmask with ones at its 2k least significant positions

            const size_t thread_count = gridDim.x * threads_per_block;
            const size_t thread_id = blockIdx.x * threads_per_block + threadIdx.x;

            /************************************************************************************************
             * Each thread is associated to a couple of reads, which it proceeds after each other.
             ***********************************************************************************************/

            for (unsigned read_id = thread_id; read_id < read_count; read_id += thread_count) {

                // construct the read from the raw sequence data
                const uint8_t* read_data_ptr = read_data + read_start_index[read_id];
                const read r(read_data_ptr, read_lengths[read_id]);

                /****************************************************************************************************
                 * Step 1: Compute the pseudo distance to *all* barcodes.
                 ***************************************************************************************************/

                // initialize the pseudo distance array
                for (unsigned barcode_id = 0; barcode_id < barcode_count; barcode_id++)
                    pseudo_distances[barcode_id * thread_count + thread_id] = 0;

                // initialize a list of barcode candidates
                unsigned barcode_candidate_count = 0;

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

                        unsigned start = k_mer_start_index[j * k_mer_count + m];
                        unsigned stop = k_mer_start_index[j * k_mer_count + m + 1];

                        // for each barcodes in which m occurs at position j
                        for (unsigned l = start; l < stop; l++) {
                            unsigned barcode_id = k_mer_map_entries[l];
                            if (pseudo_distances[barcode_id * thread_count + thread_id] == 0) {
                                // append barcode_id at the end of the threads candidate list
                                barcode_candidates[barcode_candidate_count * thread_count + thread_id] = barcode_id;
                                barcode_candidate_count++;
                            }
                            pseudo_distances[barcode_id * thread_count + thread_id] += abs(i - j) - BARCODE_LENGTH;
                        }
                    }
                }

                /*if (read_id == 1) {
                    printf("read %u has %u candidates, blockIdx.x=%u, threadIdx.x=%u\n", read_id, barcode_candidate_count, blockIdx.x, threadIdx.x);
                }*/

                /****************************************************************************************************
                 * Step 2: Determine the MAX_INDEX_SIZE barcodes with smallest pseudo distance.
                 ***************************************************************************************************/

                assert(MAX_INDEX_SIZE + 1 < barcode_count);

                struct index_element {
                    unsigned barcode_id;
                    int32_t distance;
                } my_local_index[MAX_INDEX_SIZE + 1];

                unsigned my_local_index_size = 0;

                // for each barcode id in each threads candidate list
                for (unsigned p = 0; p < barcode_candidate_count; p++) {

                    // select the p-th candidate in each threads candidate list and its associated pseudo distance
                    unsigned barcode_id = barcode_candidates[p * thread_count + thread_id];
                    int32_t q = pseudo_distances[barcode_id * thread_count + thread_id];

                    // insert {b,q} into the sorted index [0...my_index_size-1]
                    unsigned j = my_local_index_size;
                    for (; j > 0 && my_local_index[j - 1].distance > q; j--) {
                        if (j < MAX_INDEX_SIZE)
                            my_local_index[j] = my_local_index[j - 1];
                    }

                    if (j < MAX_INDEX_SIZE)
                        my_local_index[j] = {barcode_id, q};

                    if (my_local_index_size < MAX_INDEX_SIZE)
                        my_local_index_size++;
                }

                /*if (read_id == 1) {
                    printf("index of read %i: after step 2\n", read_id);
                    for (unsigned l = 0; l < my_local_index_size; l++) {
                        index_element el = my_local_index[l];
                        printf("pos=%i, barcode=%i, distance=%i\n", l, el.barcode_id, el.distance);
                    }
                }*/

                /****************************************************************************************************
                 * Step 3: Replace the pseudo distances by the SL-distance.
                 ***************************************************************************************************/

                for (unsigned rank = 0; rank < MAX_INDEX_SIZE && rank < barcode_count; rank++) {
                    barcode b = barcodes[my_local_index[rank].barcode_id];
                    my_local_index[rank].distance = distance_function::evaluate(b, r, costs);
                }

                /*if (read_id == 1) {
                    printf("index of read %i: after step 3\n", read_id);
                    for (unsigned l = 0; l < my_local_index_size; l++) {
                        index_element el = my_local_index[l];
                        printf("pos=%i, barcode=%i, distance=%i\n", l, el.barcode_id, el.distance);
                    }
                }*/

                /****************************************************************************************************
                 * Step 4: Sort the index elements by their SL-distance.
                 ***************************************************************************************************/

                // sort elements via insertion-sort
                for (unsigned i = 1; i < my_local_index_size; i++) {
                    // insert my_local_index[i] into the sorted my_local_index_size[0...i-1]
                    index_element y = my_local_index[i];
                    unsigned j = i;
                    for (; j > 0; j--) {
                        const index_element& x = my_local_index[j - 1];
                        if (x.distance > y.distance)
                            my_local_index[j] = my_local_index[j - 1];
                        else
                            break;
                    }
                    my_local_index[j] = y;
                }

                /*if (read_id == 1) {
                    printf("index of read %i: after step 4\n", read_id);
                    for (unsigned l = 0; l < my_local_index_size; l++) {
                        index_element el = my_local_index[l];
                        printf("pos=%i, barcode=%i, distance=%i\n", l, el.barcode_id, el.distance);
                    }
                }*/

                /****************************************************************************************************
                 * Step 5: Store the final index in the global memory.
                 ***************************************************************************************************/

                // write the two barcodes with smallest distance to the global memory
                out_barcodes[read_id] = my_local_index[0].barcode_id;
                out_barcodes[read_id + read_count] = my_local_index[1].barcode_id;
                out_distances[read_id] = my_local_index[0].distance;
                out_distances[read_id + read_count] = my_local_index[1].distance;
            }

        }
    }


    template <unsigned k>
    class k_mer_filtered_calling_gpu_v4 : public barcode_calling_algorithm {

    public:
        explicit k_mer_filtered_calling_gpu_v4(const distance_measure& dist, unsigned rejection_threshold)
            : barcode_calling_algorithm(std::to_string(k) + "_mer_filtered_calling_gpu_v4",
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
                  unsigned threads_per_block = 32,
                  unsigned PSEUDO_DISTANCE_WINDOW_SIZE = DEFAULT_PSEUDO_DISTANCE_WINDOW_SIZE,
                  unsigned MAX_INDEX_SIZE = DEFAULT_CANDIDATE_COUNT
        >
        extended_barcode_assignment run(
            const barcode_set& barcodes,
            const read_set& reads,
            const alignment_costs& costs) const {

            extended_barcode_assignment ass(reads.size()); // will be returned
            cudaError_t err;

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
            CUDA_CHECK(cudaMemcpy(barcodes_dev, &barcodes[0],
                barcodes.size() * sizeof(barcode), cudaMemcpyHostToDevice));

            // read start positions
            CUDA_CHECK(cudaMemcpy(read_start_index_dev, reads.get_start_indices().data(),
                reads.get_start_indices().size() * sizeof(size_t),
                cudaMemcpyHostToDevice));

            // read sequence lengths
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

            /*******************************************************************************************************************
             * Prepare the k-mer data structure
             *
             * For each combination of a k-mer x in {A,C,G,T}^k and position i in [0...sequence::LENGTH-1], we construct
             * a list of barcodes b in which x occurs in b at position i.
             *
             * Thus, k_mer_map[x][i] is a list of exactly those barcode id's b, for which x == b[i...i+k-1].
             *******************************************************************************************************************/

            k_mer_map_host_v1<k> k_mer_map(barcodes);
            unsigned k_mer_count = k_mer_map.get_k_mer_count();

            /*******************************************************************************************************************
             * Copy the k-mer data structure to GPU memory.
             *
             * For this purpose, we concatenate the individual lists of the k-mer map into a single array, similarly to an
             * adjacency array of graphs. The barcode indices associated to k_map_map[x][i] are then located at the positions
             *
             *     start_index[i * k_mer_count + x] ... start_index[i * k_mer_count + x + 1] - 1
             *
             * in this array.
             *******************************************************************************************************************/

            // k-mer data structure
            unsigned* k_mer_map_entries_dev;
            unsigned* k_mer_map_start_index_dev;

            std::vector<unsigned> k_mer_map_start_index_host(k_mer_count * (BARCODE_LENGTH - k + 1) + 1);
            std::vector<unsigned> k_mer_map_entries_host(barcodes.size() * (BARCODE_LENGTH - k + 1));
            CUDA_CHECK(cudaMalloc(&k_mer_map_entries_dev,
                barcodes.size() * (BARCODE_LENGTH - k + 1) * sizeof(unsigned)));
            CUDA_CHECK(cudaMalloc(&k_mer_map_start_index_dev,
                (k_mer_count * (BARCODE_LENGTH - k + 1) + 1) * sizeof(unsigned)));

            unsigned running_index = 0;
            for (unsigned i = 0; i <= BARCODE_LENGTH - k; i++) {

                // for each possible k-mer
                for (unsigned x = 0; x < k_mer_count; x++) {
                    // concat the lists into a single permutation of the barcodes
                    memcpy(&k_mer_map_entries_host[running_index], &k_mer_map[x][i][0],
                           k_mer_map[x][i].size() * sizeof(unsigned));
                    running_index += k_mer_map[x][i].size();
                    k_mer_map_start_index_host[i * k_mer_count + x + 1] = running_index;
                }
            }

            // copy the barcode id permutation and start indices to gpu memory
            CUDA_CHECK(cudaMemcpy(k_mer_map_entries_dev, k_mer_map_entries_host.data(),
                barcodes.size() * (BARCODE_LENGTH - k + 1) * sizeof(unsigned),
                cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(k_mer_map_start_index_dev, k_mer_map_start_index_host.data(),
                (k_mer_count * (BARCODE_LENGTH - k + 1) + 1) * sizeof(unsigned),
                cudaMemcpyHostToDevice));

            /************************************************************************************************
             * Allocate the working memory for the calling kernel.
             *
             * Each thread is responsible for a couple of reads, which it proceeds sequentially.
             *
             * The kernel computes for each read its pseudo distance to all barcodes.
             * To this end, *each* thread requires an array of size barcode_count to store the pseudo distances
             * and a second array to store candidate barcodes whose pseudo distance is not zero. As we do not
             * known how many candidate barcodes we get, this second array must be of size barcode_count, too.
             *
             * Thus, the number of CUDA threads and (as we want a constant number of threads per block) the
             * number of CUDA blocks is limited by the amount of free GPU memory.
             ************************************************************************************************/

            // working memory
            unsigned* barcode_candidates_dev;
            int32_t* pseudo_distances_dev;

            // determine the amount of free memory in bytes
            size_t free_mem = 0;
            size_t total_mem = 0;
            err = cudaMemGetInfo(&free_mem, &total_mem);
            if (err != cudaSuccess) {
                std::cerr << "Error while reading GPU memory: " << cudaGetErrorString(err) << std::endl;
                return ass;
            }

            size_t bytes_per_thread = barcodes.size() * (sizeof(unsigned) + sizeof(int32_t));
            size_t bytes_per_block = bytes_per_thread * threads_per_block;
            size_t bytes_to_use = 0.9 * free_mem; // use 90% of the free memory
            unsigned block_count = bytes_to_use / bytes_per_block;
            if (block_count * threads_per_block > reads.size())
                block_count = (reads.size() + threads_per_block - 1) / threads_per_block;
            unsigned thread_count = threads_per_block * block_count;

            CUDA_CHECK(cudaMalloc(&pseudo_distances_dev, thread_count * barcodes.size() * sizeof(int32_t)));
            CUDA_CHECK(cudaMalloc(&barcode_candidates_dev, thread_count * barcodes.size() * sizeof(unsigned)));

            /************************************************************************************************
             * Call the calling kernel.
             ************************************************************************************************/

            k_mer_filtered_calling_gpu_v4_detail::kernel
                <distance_measure, k, threads_per_block, PSEUDO_DISTANCE_WINDOW_SIZE, MAX_INDEX_SIZE>
                <<<block_count, threads_per_block>>>
                (barcodes_dev, read_start_index_dev, read_sequence_data_dev, read_lengths_dev,
                 k_mer_map_entries_dev, k_mer_map_start_index_dev,
                 pseudo_distances_dev, barcode_candidates_dev,
                 costs, barcodes.size(), reads.size(),
                 out_barcodes_dev, out_distances_dev);

            cudaDeviceSynchronize();
            err = cudaGetLastError();
            if (err != cudaSuccess) {
                std::cerr << cudaGetErrorString(err) << std::endl;
                return ass;
            }

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
            CUDA_CHECK(cudaFree(k_mer_map_start_index_dev));
            CUDA_CHECK(cudaFree(pseudo_distances_dev));
            CUDA_CHECK(cudaFree(barcode_candidates_dev));

            return ass;
        }


    };


    ;


}


#endif //INC_2OPT_KMER_FILTERED_CALLING_GPU_V4_CUH
