//
// Created by steffen on 25.06.24.
//

#ifndef INC_2OPT_KMER_FILTERED_CALLING_GPU_V4_CUH
#define INC_2OPT_KMER_FILTERED_CALLING_GPU_V4_CUH

#include "barcode.h"
#include "read.h"
#include "barcode_assignment.h"
#include <vector>


template <
    unsigned k, // length of k-mers
    unsigned threads_per_block,
    uint8_t (*exact_distance)(const sequence&, const sequence&, uint8_t upper_bound) // distance measure
>
__global__ void calling_kernel_kmer_filter_gpu_v4(
        const barcode *barcodes,
        const read *reads,
        const unsigned *k_mer_map_entries,  // k-mer data structure
        const unsigned *k_mer_start_index,  // k-mer data structure
        unsigned *index_barcodes,           // for each block: array of size MAX_INDEX_SIZE
        uint8_t *index_distances,           // for each block: array of size MAX_INDEX_SIZE
        int16_t *pseudo_distances,          // for each block: array of size barcode_count
        unsigned *barcode_candidates,
        const unsigned barcode_count,       // number of barcodes
        const unsigned read_count)          // number of reads
{
    static constexpr uint32_t k_mer_count = 1 << 2 * k; // 4^k
    static constexpr uint32_t mask = (1 << 2 * k) - 1; // bitmask with ones at its 2k least significant positions

    const size_t thread_count = gridDim.x * threads_per_block;
    const size_t thread_id = blockIdx.x * threads_per_block + threadIdx.x;

    /************************************************************************************************
     * Each thread is associated to a couple of reads, which it proceeds after each other.
     ***********************************************************************************************/

    for(unsigned read_id = thread_id; read_id < read_count; read_id += thread_count) {

        const read r = reads[read_id];

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
            m = (m << 2) + r[l];

        // for each position i in which a k-mer may start
        for (int i = 0; i + k <= sequence::LENGTH; i++) {

            // invariant: m == r[i...i+k-1]
            m = ((m << 2) & mask) + r[i + k - 1];

            // for each position j in a window around position i
            int j = 0;
            if (i > PSEUDO_DISTANCE_WINDOW_SIZE)
                j = i - PSEUDO_DISTANCE_WINDOW_SIZE;
            for (; j + k <= sequence::LENGTH && j <= i + PSEUDO_DISTANCE_WINDOW_SIZE; j++) {

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
                    pseudo_distances[barcode_id * thread_count + thread_id] += abs(i - j) - sequence::LENGTH;
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
            int16_t distance;
        } my_local_index[MAX_INDEX_SIZE + 1];

        unsigned my_local_index_size = 0;

        // for each barcode id in each threads candidate list
        for (unsigned p = 0; p < barcode_candidate_count; p++) {

            // select the p-th candidate in each threads candidate list and its associated pseudo distance
            unsigned barcode_id = barcode_candidates[p * thread_count + thread_id];
            int16_t q = pseudo_distances[barcode_id * thread_count + thread_id];

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
            my_local_index[rank].distance = exact_distance(b, r, UINT8_MAX);
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
                const index_element &x = my_local_index[j - 1];
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

        for (unsigned rank = 0; rank < MAX_INDEX_SIZE && rank < barcode_count; rank++) {
            index_barcodes[rank * read_count + read_id] = my_local_index[rank].barcode_id;
            index_distances[rank * read_count + read_id] = my_local_index[rank].distance;
        }
    }
}




struct k_mer_filtered_calling_gpu_v4 : public barcode_assignment {

    k_mer_filtered_calling_gpu_v4(const std::vector<barcode> &barcodes,
                                  const std::vector<read> &reads,
                                  const unsigned k,
                                  const int distance_measure,
                                  const int rejection_threshold = REJECTION_THRESHOLD)
            : barcode_assignment(run(barcodes, reads, k, distance_measure, rejection_threshold)) {}

private:

    barcode_assignment run(const std::vector<barcode>& barcodes,
                                  const std::vector<read>& reads,
                                  const unsigned k,
                                  const int distance_measure,
                                  const int rejection_threshold) {
        switch (k) {
        case 4:
            return run<4>(barcodes, reads, distance_measure, rejection_threshold);
        case 5:
            return run<5>(barcodes, reads, distance_measure, rejection_threshold);
        case 6:
            return run<6>(barcodes, reads, distance_measure, rejection_threshold);
        case 7:
            return run<7>(barcodes, reads, distance_measure, rejection_threshold);
        default:
            std::cerr << "invalid parameter k!" << std::endl;
            return barcode_assignment(reads.size()); // return invalid assignment
        }
    }

    template <unsigned k>
    barcode_assignment run(const std::vector<barcode>& barcodes,
                              const std::vector<read>& reads,
                              const int distance_measure,
                              const int rejection_threshold) {
        switch (distance_measure) {
        case SEQUENCE_LEVENSHTEIN_DISTANCE:
            return run<k, sequence_levenshtein_distance_v5>(barcodes, reads, rejection_threshold);
        case LEVENSHTEIN_DISTANCE:
            return run<k, levenshtein_distance_v2>(barcodes, reads, rejection_threshold);
        default:
            std::cerr << "invalid distance measure!" << std::endl;
            return barcode_assignment(reads.size()); // return invalid assignment
        }
    }

    template <
        unsigned k,
        uint8_t (*exact_distance)(const sequence&, const sequence&, uint8_t upper_bound) // distance measure
    >
    barcode_assignment run(const std::vector<barcode> &barcodes,
                                const std::vector<read> &reads,
                                const int rejection_threshold) {

        barcode_assignment ass(reads.size()); // will be returned

        /************************************************************************************************
         * Allocate gpu memory for the barcodes and reads and copy them to the GPU.
         ************************************************************************************************/

        barcode *barcodes_dev;              // barcode array on gpu memory
        read *reads_dev;                    // read array on gpu memory

        cudaMalloc(&barcodes_dev, barcodes.size() * sizeof(barcode));
        cudaMalloc(&reads_dev, reads.size() * sizeof(read));
        cudaMemcpy(barcodes_dev, &barcodes[0], barcodes.size() * sizeof(barcode), cudaMemcpyHostToDevice);
        cudaMemcpy(reads_dev, &reads[0], reads.size() * sizeof(read), cudaMemcpyHostToDevice);

        /************************************************************************************************
         * Allocate the output memory for the calling kernel.
         *
         * For each read r, the kernel produces a list of at most MAX_INDEX_SIZE tuples of the form
         *
         *    (barcode_id b, pseudo_distance q(r,b), SL-dist(r,b)).
         ************************************************************************************************/

        unsigned *out_barcodes_dev;
        int16_t *out_pseudo_distances_dev;
        uint8_t *out_sl_distances_dev;
        unsigned *barcode_candidates_dev;

        cudaMalloc(&out_barcodes_dev, reads.size() * MAX_INDEX_SIZE * sizeof(unsigned)); // todo: cudaMallocPitch
        cudaMalloc(&out_sl_distances_dev, reads.size() * MAX_INDEX_SIZE * sizeof(uint8_t)); // todo: cudaMallocPitch

        /*******************************************************************************************************************
         * For each combination of a k-mer x in {A,C,G,T}^k and position i in [0...sequence::LENGTH-1], we construct
         * a list of barcodes b in which x occurs in b at position i.
         *
         * Thus, k_mer_map[x][i] is a list of exactly those barcode id's b, for which x == b[i...i+k-1].
         *******************************************************************************************************************/

        static_assert(2 * k < 32, "invalid value of k");
        static constexpr uint32_t k_mer_count = 1 << 2 * k; // 4^k

        // create empty k-mer-map
        std::vector<std::vector<std::vector<unsigned>>> k_mer_map(
                k_mer_count, std::vector<std::vector<unsigned >>(sequence::LENGTH));

        // for each k-mer position i
#pragma omp parallel for
        for (unsigned i = 0; i <= sequence::LENGTH - k; i++) {

            // for each barcode
            for (unsigned barcode_id = 0; barcode_id < barcodes.size(); barcode_id++) {
                const barcode &b = barcodes[barcode_id];

                // consider the k-mer x = b[i...i+k-1]
                uint32_t x = 0; // x = b[i...i+k-1]
                for (unsigned l = 0; l < k; l++)
                    x = (x << 2) + b[i + l];
                assert(x < k_mer_count);

                // k-mer x is located in barcode b at position i
                k_mer_map[x][i].push_back(barcode_id);
            }
        }

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

        unsigned *map_entries_dev, *map_entries_host;
        unsigned *start_index_dev, *start_index_host;

        cudaMalloc(&map_entries_dev, barcodes.size() * (sequence::LENGTH - k + 1) * sizeof(unsigned));
        cudaMalloc(&start_index_dev, (k_mer_count * (sequence::LENGTH - k + 1) + 1) * sizeof(unsigned));
        start_index_host = new unsigned[k_mer_count * (sequence::LENGTH - k + 1) + 1];
        map_entries_host = new unsigned[barcodes.size() * (sequence::LENGTH - k + 1)];

        unsigned running_index = 0;
        for (unsigned i = 0; i <= sequence::LENGTH - k; i++) {

            // for each possible k-mer
            for (unsigned x = 0; x < k_mer_count; x++) {
                // concat the lists into a single permutation of the barcodes
                memcpy(&map_entries_host[running_index], &k_mer_map[x][i][0],
                       k_mer_map[x][i].size() * sizeof(unsigned));
                running_index += k_mer_map[x][i].size();
                start_index_host[i * k_mer_count + x + 1] = running_index;
            }
        }

        // copy the barcode id permutation and start indices to gpu memory
        cudaMemcpy(map_entries_dev, map_entries_host,
                   barcodes.size() * (sequence::LENGTH - k + 1) * sizeof(unsigned),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(start_index_dev, start_index_host,
                   (k_mer_count * (sequence::LENGTH - k + 1) + 1) * sizeof(unsigned),
                   cudaMemcpyHostToDevice);
        delete[] start_index_host;
        delete[] map_entries_host;

        /********************************************************************************************************
         * As each CUDA thread needs two arrays of size barcode_count each, the number of CUDA blocks is limited
         * by the amount of free GPU memory.
         ******************************************************************************************************/

        size_t free_mem = 0;
        size_t total_mem = 0;
        cudaError_t cuda_status = cudaMemGetInfo(&free_mem, &total_mem);
        if (cuda_status != cudaSuccess) {
            std::cerr << "Error while reading GPU memory: " << cudaGetErrorString(cuda_status) << std::endl;
            // todo: free GPU memory
            return barcode_assignment(reads.size());
        }

        static constexpr unsigned threads_per_block = 32;
        unsigned block_count = (unsigned) (0.9 * free_mem /
                                           ((double) threads_per_block * barcodes.size() * (sizeof(int16_t) + sizeof(unsigned))));
        if (block_count * threads_per_block > reads.size())
            block_count = (reads.size() + threads_per_block - 1) / threads_per_block;
        unsigned thread_count = threads_per_block * block_count;

        //std::cout << "block_count = " << block_count << std::endl;

        cudaMalloc(&out_pseudo_distances_dev,
                   thread_count * barcodes.size() * sizeof(int16_t)); // todo: cudaMallocPitch
        cudaMalloc(&barcode_candidates_dev,
                   thread_count * barcodes.size() * sizeof(unsigned)); // todo: cudaMallocPitch

        /************************************************************************************************
         * Call the calling kernel.
         ************************************************************************************************/

        calling_kernel_kmer_filter_gpu_v4
        <k, threads_per_block, exact_distance>
        <<<block_count, threads_per_block>>>(
                barcodes_dev, reads_dev, map_entries_dev, start_index_dev, out_barcodes_dev,
                        out_sl_distances_dev, out_pseudo_distances_dev, barcode_candidates_dev, barcodes.size(), reads.size()
        );

        auto err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << cudaGetErrorString(err) << std::endl;
            return barcode_assignment(reads.size());
        }

        /************************************************************************************************
         * Copy the barcode indices to the host memory and select the first of each barcode.
         ************************************************************************************************/

        uint8_t *out_sl_distances_host = new uint8_t[reads.size() * MAX_INDEX_SIZE];

        cudaMemcpy(&ass[0], out_barcodes_dev, reads.size()*sizeof(unsigned), cudaMemcpyDeviceToHost );
        cudaMemcpy(out_sl_distances_host, out_sl_distances_dev, reads.size() * sizeof(uint8_t), cudaMemcpyDeviceToHost);

#pragma omp parallel for
        for(unsigned read_id = 0; read_id < reads.size(); read_id++){
            if(out_sl_distances_host[read_id] > rejection_threshold)
                ass[read_id] = UINT_MAX;
        }

        /************************************************************************************************
         * Free the memory.
         ************************************************************************************************/

        cudaFree(barcodes_dev);
        cudaFree(reads_dev);
        cudaFree(out_barcodes_dev);
        cudaFree(out_sl_distances_dev);
        cudaFree(out_pseudo_distances_dev);
        cudaFree(map_entries_dev);
        cudaFree(start_index_dev);
        cudaFree(barcode_candidates_dev);
        delete[] out_sl_distances_host;

        return ass;
    }


};

#endif //INC_2OPT_KMER_FILTERED_CALLING_GPU_V4_CUH
