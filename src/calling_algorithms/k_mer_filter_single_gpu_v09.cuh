//
// Created by steffen on 12.03.26.
//

#pragma once

#include "barcode_calling_algorithm.h"
#include "../k_mer_map_dev_chunked_v2.h"
#include "../read_set_dev.cuh"
#include "../distance/weighted_levenshtein_v1.cuh"
#include "../distance/weighted_sequence_levenshtein_v1.cuh"
#include "../cuda_helper.cuh"
#include "../arg_min.cuh"
#include "../memory_alignment_helper.cuh"

namespace barcode_calling {
    namespace k_mer_filter_gpu_v9_detail {


        /**
         * Use binary search to find a position i in [start_index,CANDIDATE_COUNT-1] such that
         *
         *   a) the distances of all candidates left of i are <= dist, i.e.
         *
         *      candidate_distance[j] <= dist for all j in [start_index,i-1], and
         *
         *   b) the distances of all candidates right of i are > dist, i.e.
         *
         *      candidate_distance[j] > dist for all j in [i,CANDIDATE_COUNT-1].
         *
         * The running time of this function is log2(CANDIDATE_COUNT-start_index).
         *
         * @param dist New distance to be inserted.
         * @param candidate_distance CANDIDATE_COUNT (pseudo) distances in ascending order.
         * @param start_index
         * @return
         */
        template <unsigned CANDIDATE_COUNT>
        __device__ unsigned find_insert_position(
            const int32_t dist,
            const int32_t* candidate_distance,
            const unsigned start_index = 0) {

            unsigned left = start_index;
            unsigned right = CANDIDATE_COUNT;

            /***********************************************************************************
             * Invariant:
             *
             *  a) for all j in [start_index,left-1]:      candidate_distance[j] <= dist
             *  b) for all j in [right,CANDIDATE_COUNT-1]: candidate_distance[j] >  dist
             *
             * Both invariantes are trivially fulfilled when entering the loop.
             *
             * When left == right:
             *
             * For all j in [start_index,left-1]:     candidate_distance[j] <= dist
             * For all j in [left,CANDIDATE_COUNT-1]: candidate_distance[j] >  dist
             *
             * Thus, we can return left as the result.
             ***********************************************************************************/

            for (unsigned mid = (left + right) / 2; left < right; mid = (left + right) / 2) {
                assert(mid < CANDIDATE_COUNT);
                assert(mid >= start_index);
                if (candidate_distance[mid] <= dist)
                    left = mid + 1;
                else // candidate_distance[mid] > dist
                    right = mid;
            }
            return left;
        }


        template <typename distance_function,
                  unsigned k, // length of k-mers
                  unsigned barcodes_per_chunk,
                  unsigned threads_per_block,
                  unsigned PSEUDO_DISTANCE_WINDOW_SIZE,
                  unsigned CANDIDATE_COUNT
        >
        __global__ void kernel(

            // input variables
            const barcode* barcodes, // linear array of all barcodes
            const size_t barcode_count, // total number of all barcodes
            const read_set_dev reads,
            const k_mer_map_dev_chunked_v2<k, barcodes_per_chunk> k_mer_maps, // for each barcode chunk one map
            const distance_function distance,

            // output variables
            unsigned* out_barcodes, // out: 2*read_count entries
            int32_t* out_distances // out: 2*read_count entries
        ) {

            //auto start = clock64();

            const unsigned chunk_count = (barcode_count + barcodes_per_chunk - 1) / barcodes_per_chunk;

            assert(CANDIDATE_COUNT <= barcode_count);
            assert(barcodes_per_chunk >= CANDIDATE_COUNT);

            /******************************************************************************************
             * Each block is associated to a single read r.
             *
             * In the first step, the block determines a list of barcodes with smallest pseudo distance.
             * These barcodes are candidates for the final assignment.
             *****************************************************************************************/

            const unsigned read_id = blockIdx.x;
            assert(read_id < reads.size());
            const read r = reads[read_id];

            assert(r.length() >= k);

            /*if (threadIdx.x == 0 && read_id < 10) {

                char str[200];
                r.bytes_to_string(str);
                printf("thread_id=%u, read_id=%u, read_len=%u, read=%s, r.sequence_data=%p, *r.sequence_data=%u\n",
                    threadIdx.x, read_id, reads.get_read_length(read_id), str, r.sequence_data, *r.sequence_data);
            }
            __syncthreads();

            return;*/

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
             *       aux_barcode_ids = unsigned[max(CANDIDATE_COUNT,threads_per_block)]
             *
             * 4. We use an auxiliary array to temporarily store distances.
             *
             *       aux_distances = int32_t[CANDIDATE_COUNT]
             *********************************************************************************/

            extern __shared__ unsigned char smem[];

            memory_alignment_helper aligner(smem);
            auto pseudo_distance = aligner.register_array<int32_t>(barcodes_per_chunk);
            auto candidate_barcodes = aligner.register_array<unsigned>(CANDIDATE_COUNT + 1);
            auto candidate_distance = aligner.register_array<int32_t>(CANDIDATE_COUNT + 1);
            auto aux_barcode_ids = aligner.register_array<unsigned>(
                max(CANDIDATE_COUNT, threads_per_block));
            auto aux_distances = aligner.register_array<int32_t>(CANDIDATE_COUNT);

            /****************************************************************************
             * initialize the candidate list with dummy elements
             ***************************************************************************/

            for (unsigned i = threadIdx.x; i <= CANDIDATE_COUNT; i += threads_per_block) {
                candidate_barcodes[i] = UINT_MAX;
                candidate_distance[i] = INT32_MAX;
            }
            __syncthreads();

            /****************************************************************************
             * For each chunk of barcodes, one after the other
             ***************************************************************************/

            for (unsigned chunk_id = 0; chunk_id < chunk_count; ++chunk_id) {

                /****************************************************************************
                 * Step 1: Compute the pseudo distance to all barcodes of the current chunk.
                 ***************************************************************************/

                // initialize the pseudo distance array
                for (unsigned barcode_id = threadIdx.x; barcode_id < barcodes_per_chunk;
                    barcode_id += threads_per_block) {
                    pseudo_distance[barcode_id] = 0;
                }

                __syncthreads();

                const uint32_t mask = (1 << 2 * k) - 1; // bitmask with ones at its 2k least significant positions

                // prepare the first k-1 characters of k-mer m = r[0...k-1]
                uint32_t m = 0; // m = r[0...k-2]
                for (unsigned l = 0; l < k - 1; l++)
                    m = (m << 2) + r[l].to_uint8();

                // for each position i in which a k-mer may start in the read
                for (int i = 0; i + k <= BARCODE_LENGTH && i + k <= r.length(); i++) {

                    // invariant: m == r[i...i+k-1]
                    m = ((m << 2) & mask) + r[i + k - 1].to_uint8();

                    /***********************************************************************
                     * Consider a window of at most 2*PSEUDO_DISTANCE_WINDOW_SIZE+1
                     * positions around i. The window is given by [left, right].
                     ***********************************************************************/

                    int left = 0;
                    if (i > PSEUDO_DISTANCE_WINDOW_SIZE)
                        left = i - PSEUDO_DISTANCE_WINDOW_SIZE;
                    int right = i + PSEUDO_DISTANCE_WINDOW_SIZE;
                    if (right + k > BARCODE_LENGTH)
                        right = BARCODE_LENGTH - k;

                    /***********************************************************************
                     * We determine the list of all barcodes in which m starts some position
                     * within [left,right].
                     *
                     * For each barcode in this list, we get a second list of the same length
                     * in which the starting positions of the associated k-mer m are stored.
                     ***********************************************************************/

                    auto barcodes_begin = k_mer_maps.get_barcode_ids(chunk_id, m, left);
                    auto barcodes_end = k_mer_maps.get_barcode_ids(chunk_id, m, right + 1);
                    auto k_mer_positions_begin = k_mer_maps.get_k_mer_start_positions(chunk_id, m, left);
                    auto k_mer_positions_end = k_mer_maps.get_k_mer_start_positions(chunk_id, m, right + 1);

                    size_t barcode_list_length = barcodes_end - barcodes_begin;
                    assert(barcode_list_length == k_mer_positions_end - k_mer_positions_begin);

                    // for each barcode in which m starts at some position within [left,right]
                    for (size_t l = threadIdx.x; l < barcode_list_length; l += threads_per_block) {

                        unsigned global_barcode_id = *(barcodes_begin + l);
                        unsigned local_barcode_id = global_barcode_id % barcodes_per_chunk;
                        uint8_t j = *(k_mer_positions_begin + l);
                        //printf("thread=%i, chunk_id=%u, i=%i, m=%u, left=%u, right=%u, barcode_list_length=%lu, l=%lu, barcode_id=%u, j=%u\n", threadIdx.x, chunk_id, i, m, left, right, barcode_list_length, l, barcode_id, j);

                        assert(j >= left);
                        assert(j <= right);

                        /******************************************************************************
                         * The k-mer m starts at position j in the barcode with the given id.
                         *
                         * Increase the pseudo distance of this barcode by |i-j| - BARCODE_LENGTH.
                         ******************************************************************************/

                        // atomic operation needed as barcodes may occur twice in the same list
                        atomicAdd(&pseudo_distance[local_barcode_id], abs(i - j) - BARCODE_LENGTH);
                    }

                    __syncthreads();
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
                     ****************************************************************************************/

                    unsigned my_best_barcode_id = threadIdx.x;
                    for (unsigned barcode_id = threadIdx.x + threads_per_block;
                         barcode_id < barcodes_per_chunk; barcode_id += threads_per_block) {
                        if (pseudo_distance[my_best_barcode_id] > pseudo_distance[barcode_id])
                            my_best_barcode_id = barcode_id;
                    }

                    /*****************************************************************************************
                     * Find a barcode with minimum pseudo distance within the current chunk.
                     ****************************************************************************************/

                    aux_barcode_ids[threadIdx.x] = my_best_barcode_id;
                    __threadfence_block();
                    __syncthreads();

                    unsigned best_barcode_in_chunk = arg_min<threads_per_block>(
                        aux_barcode_ids, pseudo_distance, threads_per_block);

                    //printf("thread=%i, my_best_barcode_id=%u, best_barcode_in_chunk=%u\n", threadIdx.x, my_best_barcode_id, best_barcode_in_chunk);
                    assert(best_barcode_in_chunk < barcodes_per_chunk);
                    int32_t best_pseudo_dist_in_chunk = pseudo_distance[best_barcode_in_chunk];

                    __syncthreads();

                    /*****************************************************************************************
                     * NEW: Use binary search to find the position to which the best barcode in this chunk
                     * must be inserted.
                     ****************************************************************************************/

                    const unsigned insert_position = find_insert_position<CANDIDATE_COUNT>
                        (best_pseudo_dist_in_chunk, candidate_distance, rank);
                    assert(insert_position <= CANDIDATE_COUNT);
                    assert(insert_position >= rank);

                    /*****************************************************************************************
                     * We may break the loop if the best barcode in the current chunk is worst than all
                     * candidate barcodes.
                     ****************************************************************************************/

                    if (insert_position == CANDIDATE_COUNT)
                        break;

                    /*****************************************************************************************
                     * NEW: All threads collectively move the candidates right of the insertion position
                     * to the auxiliary space within the shared memory.
                     ****************************************************************************************/

                    for (unsigned i = insert_position + threadIdx.x; i + 1 < CANDIDATE_COUNT; i += threads_per_block) {
                        aux_barcode_ids[i] = candidate_barcodes[i];
                        aux_distances[i] = candidate_distance[i];
                    }
                    __syncthreads();

                    /*****************************************************************************************
                     * Thread 0 inserts the best barcode in the current chunk to the right candidate position.
                     *
                     * It sets its pseudo distance to infinity so it cannot be used a second time.
                     ****************************************************************************************/

                    if (threadIdx.x == 0) {
                        candidate_barcodes[insert_position] = chunk_id * barcodes_per_chunk + best_barcode_in_chunk;
                        candidate_distance[insert_position] = best_pseudo_dist_in_chunk;
                        pseudo_distance[best_barcode_in_chunk] = INT32_MAX;
                    }
                    __syncthreads();

                    /*****************************************************************************************
                     * NEW: All threads collectively move the candidate barcodes and distances back to their
                     * new right place, which is one position right of their previous one.
                     ****************************************************************************************/

                    for (unsigned i = insert_position + threadIdx.x; i + 1 < CANDIDATE_COUNT; i += threads_per_block) {
                        candidate_barcodes[i + 1] = aux_barcode_ids[i];
                        candidate_distance[i + 1] = aux_distances[i];
                    }
                    __syncthreads();
                    __threadfence_block();
                }
            }

            /*auto phase1 = clock64() - start;
            start = clock64();*/

            /****************************************************************************************************
             * Step 3: Replace the pseudo distances by the exact distance.
             ***************************************************************************************************/

            for (unsigned i = threadIdx.x; i < CANDIDATE_COUNT && candidate_barcodes[i] < barcode_count;
                 i += threads_per_block) {
                barcode b = barcodes[candidate_barcodes[i]];
                candidate_distance[i] = distance(b, r);
            }

            __syncthreads();

            /* auto phase2 = clock64() - start;
             start = clock64();*/

            /*********************************************************************************
             * Step 4: Determine the barcode with minimum exact distance.
             *********************************************************************************/

            for (unsigned i = threadIdx.x; i < CANDIDATE_COUNT; i += threads_per_block)
                aux_barcode_ids[i] = i;

            __syncthreads();

            // in-block reduction
            unsigned i = arg_min<32>(aux_barcode_ids, candidate_distance, CANDIDATE_COUNT);

            if (threadIdx.x == 0) {
                unsigned barcode_id_1 = candidate_barcodes[i];
                int32_t distance_1 = candidate_distance[i];

                // write to global memory
                out_barcodes[read_id] = barcode_id_1;
                out_distances[read_id] = distance_1;

                /* printf("device_id=%i, read_id=%u, best_barcode=%u, best_distance=%i\n",
                        device_id, read_id, barcode_id_1, distance_1);*/

                // make barcode_id_1 invalid
                candidate_distance[i] = INT32_MAX;
            }

            __syncthreads();
            __threadfence_block();

            /*********************************************************************************
             * Step 5: Determine the barcode with the second to minimum exact distance.
             *********************************************************************************/

            for (unsigned i = threadIdx.x; i < CANDIDATE_COUNT; i += threads_per_block)
                aux_barcode_ids[i] = i;
            __syncthreads();

            // in-block reduction
            unsigned j = arg_min<32>(aux_barcode_ids, candidate_distance, CANDIDATE_COUNT);

            if (threadIdx.x == 0) {
                unsigned barcode_id_2 = candidate_barcodes[j];
                int32_t distance_2 = candidate_distance[j];

                // write to global memory
                out_barcodes[read_id + reads.size()] = barcode_id_2;
                out_distances[read_id + reads.size()] = distance_2;
            }

            /*auto phase3 = clock64() - start;
            start = clock64();

            if (threadIdx.x == 0 && read_id == 0) {
                printf("phase 1: %lu\nphase 2: %lu\nphase 3: %lu\n", phase1, phase2, phase3);
            }*/

        }
    }


    template <unsigned k,
              typename distance_function,
              unsigned barcodes_per_chunk = 2048,
              unsigned threads_per_block = 64,
              unsigned PSEUDO_DISTANCE_WINDOW_SIZE = DEFAULT_PSEUDO_DISTANCE_WINDOW_SIZE,
              unsigned CANDIDATE_COUNT = DEFAULT_CANDIDATE_COUNT>
    class k_mer_filter_single_gpu_v09 : public barcode_calling_algorithm {

    protected:

        cuda_device dev;

        // gpu input data
        barcode* barcodes_dev = nullptr;
        read_set_dev reads_dev;

        // gpu output data
        unsigned* out_barcodes_dev = nullptr;
        int32_t* out_distances_dev = nullptr;

        // k-mer datastructures in device memory
        k_mer_map_dev_chunked_v2<k, barcodes_per_chunk> k_mer_map_dev;

        const distance_function distance;

    public:

        static std::string name() {
            return std::to_string(k) + "_mer_filter_single_gpu_v9<"
                + std::to_string(barcodes_per_chunk)
                + "," + std::to_string(threads_per_block) + ">";
        }

        explicit k_mer_filter_single_gpu_v09(
            const barcode_set& barcodes,
            const read_set& reads,
            const int32_t rejection_threshold = INT32_MAX,
            const distance_function dist = distance_function(),
            const cuda_device& dev = cuda_device())
            : barcode_calling_algorithm(barcodes, reads, rejection_threshold), distance(dist), dev(dev) {}

        /**
         *  Do the calculation.
         **/
        k_mer_filter_single_gpu_v09& run() override {
            return start().finalize();
        }

        /**
        * Start the algorithm asynchronously.
        */
        virtual k_mer_filter_single_gpu_v09& start() {

             /*std::cout << "k_mer_filter_single_gpu_v9::start() \n"
                 << " barcodes.size() = " << barcodes.size() << "\n"
                 << " reads.size()    = " << reads.size() << "\n"
                 << " dev.id          = " << dev.id << "\n"
                 << " dev.stream      = " << dev.stream
                 << std::endl;*/

            CUDA_CHECK(cudaSetDevice(dev.id));

            /*******************************************************************************
             * The main idea of this algorithm is to store the pseudo distances
             * within the GPU's shared memory instead of the global memory.
             * In doing so, the calculation of the pseudo distances may be faster.
             *
             * For this purpose, we divide the barcode set into chunks of
             * approximately the same size. The chunks are processed sequentially.
             *
             * The shared memory of each block is organized as described within the kernel.
             ********************************************************************************/

            memory_alignment_helper aligner(nullptr);
            aligner.register_array<int32_t>(barcodes_per_chunk);
            aligner.register_array<unsigned>(CANDIDATE_COUNT + 1);
            aligner.register_array<int32_t>(CANDIDATE_COUNT + 1);
            aligner.register_array<unsigned>(max(CANDIDATE_COUNT, threads_per_block));
            aligner.register_array<int32_t>(CANDIDATE_COUNT);
            size_t shared_size = aligner.size();

            // determine the size of the shared memory and the number of barcodes chunks
            cudaDeviceProp prop{};
            CUDA_CHECK(cudaGetDeviceProperties(&prop, dev.id));

            if (shared_size > prop.sharedMemPerBlock) {
                /*std::cout << "required shared memory per block:  " << shared_size << std::endl;
                std::cout << "available shared memory per block: " << prop.sharedMemPerBlock << std::endl;*/

                throw std::runtime_error(
                    "Insufficient shared memory on target GPU for " + std::to_string(barcodes_per_chunk) +
                    " barcodes per chunk!");
            }

            /***********************************************************************
             * Allocate gpu memory for the barcodes and move the data to the GPU
             ************************************************************************/

            CUDA_CHECK(cudaMalloc(&barcodes_dev, barcodes.size() * sizeof(barcode)));
            CUDA_CHECK(cudaMemcpyAsync(barcodes_dev, barcodes.data(),
                barcodes.size() * sizeof(barcode), cudaMemcpyHostToDevice, dev.stream));

            /************************************************************************************************
            * Copy the reads to gpu memory.
            ************************************************************************************************/

            reads_dev.init(reads, dev);

            /************************************************************************************************
             * Allocate the output memory for the calling kernel.
             *
             * For each read r, the kernel outputs the 2 closest barcodes together with the distances.
             ************************************************************************************************/

            CUDA_CHECK(cudaMalloc(&out_barcodes_dev, 2*reads.size() * sizeof(unsigned)));
            CUDA_CHECK(cudaMalloc(&out_distances_dev, 2*reads.size() * sizeof(int32_t)));

            /****************************************************************************
             * For each barcode chunk, we construct an associated k-mer map.
             *
             * For each combination of a k-mer x in {A,C,G,T}^k and position i
             * in [0...BARCODE_LENGTH-1], this map contains a list of barcodes b in which
             * x occurs in b at position i.
             *
             * Thus, k_mer_map[x][i] is a list of barcodes b for which x == b[i...i+k-1].
             ****************************************************************************/

            k_mer_map_dev.init(barcodes, dev);

            /************************************************************************************************
             * Select and call the kernel.
             ************************************************************************************************/

            k_mer_filter_gpu_v9_detail::kernel
                <distance_function, k, barcodes_per_chunk, threads_per_block,
                 PSEUDO_DISTANCE_WINDOW_SIZE, CANDIDATE_COUNT>
                <<<reads.size(), threads_per_block, shared_size, dev.stream>>>
                (barcodes_dev, barcodes.size(), reads_dev, k_mer_map_dev,
                 distance, out_barcodes_dev, out_distances_dev);

            CUDA_CHECK(cudaGetLastError());

            /****************************************************************************
             * Copy the closest barcodes indices and distances to the host memory.
             ****************************************************************************/

            CUDA_CHECK(cudaMemcpyAsync(closest_barcodes,
                out_barcodes_dev,
                reads.size() * sizeof(unsigned),
                cudaMemcpyDeviceToHost, dev.stream));

            CUDA_CHECK(cudaMemcpyAsync(closest_barcodes_2nd,
                out_barcodes_dev + reads.size(),
                reads.size() * sizeof(unsigned),
                cudaMemcpyDeviceToHost, dev.stream));

            CUDA_CHECK(cudaMemcpyAsync(closest_distances,
                out_distances_dev,
                reads.size() * sizeof(int32_t), cudaMemcpyDeviceToHost,
                dev.stream));

            CUDA_CHECK(cudaMemcpyAsync(closest_distances_2nd,
                out_distances_dev + reads.size(),
                reads.size() * sizeof(int32_t), cudaMemcpyDeviceToHost,
                dev.stream));

            return *this;
        }

        /**
         * Wait for the computation to complete.
         * @return
         */
        k_mer_filter_single_gpu_v09& finalize() {

            /****************************************************************************
            * Wait for the stream to complete its tasks.
            ****************************************************************************/

            CUDA_CHECK(cudaStreamSynchronize(dev.stream));

            /************************************************************************************************
             * Free the memory.
             ************************************************************************************************/

            CUDA_CHECK(cudaFree(barcodes_dev));
            reads_dev.finalize();
            CUDA_CHECK(cudaFree(out_barcodes_dev));
            CUDA_CHECK(cudaFree(out_distances_dev));
            k_mer_map_dev.finalize();

            return *this;
        }
    };
}
