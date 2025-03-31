//
// Created by steffen on 25.06.24.
//

#ifndef INC_2OPT_KMER_FILTERED_CALLING_GPU_V1_CUH
#define INC_2OPT_KMER_FILTERED_CALLING_GPU_V1_CUH

#include "barcode.h"
#include "read.h"
#include "barcode_assignment.h"
#include "pseudo_distance.h"
#include "index_element.h"
#include "sequence_levenshtein_distance_v2.cuh"
#include "levenshtein_distance_v2.cuh"
#include <vector>

#include "sequence_levenshtein_distance_v5.cuh"

template <
    unsigned k, // length of k-mers
    uint8_t (*exact_distance)(const sequence&, const sequence&, uint8_t upper_bound) // distance measure
>
__global__ void calling_kernel_kmer_filter_gpu_v1(
    const barcode* barcodes, // pointer to all barcodes
    const read* reads, // pointer to all reads
    index_element* index, //
    unsigned barcode_count, // number of barcodes
    unsigned read_count // number of reads
) {
    /*****************************************************************************************************
     * Each thread searches the best matching barcode for a single read r.
     ****************************************************************************************************/

    unsigned read_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (read_id >= read_count)
        return;

    read r = reads[read_id];

    /*****************************************************************************************************
     * To find a best matching barcode for the read r, the thread first determines a set of MAX_INDEX_SIZE
     * barcodes that have minimum pseudo distance to r. This set of barcodes is called index.
     *
     * Each index position stores a pair {b,q} of barcode id and its associated pseudo distance.
     * The index entries are sorted increasingly by their pseudo distance.
     ****************************************************************************************************/

    unsigned my_index_size = 0; // number of elements in the index

    // for each barcode b
    for (unsigned barcode_id = 0; barcode_id < barcode_count; barcode_id++) {
        barcode b = barcodes[barcode_id];

        // calculate pseudo distance between b and r
        int16_t q = pseudo_distance<k, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);

        // insert {b,q} into the sorted index [0...my_index_size-1]

        // probe element j-1 of the index associated to read r
        unsigned j = my_index_size;
        for (; j > 0 && index[(j - 1) * read_count + read_id].distance > q; j--) {
            if (j < MAX_INDEX_SIZE)
                index[j * read_count + read_id] = index[(j - 1) * read_count + read_id];
        }
        if (j < MAX_INDEX_SIZE)
            index[j * read_count + read_id] = {barcode_id, q};

        if (my_index_size < MAX_INDEX_SIZE)
            my_index_size++;
    }

    /**************************************************************************************************************
     * Now the index contains the barcodes with smallest pseudo distance to the read r.
     *
     * Next, we replace each pseudo-distance by the sequence Levenshtein distance (SL-distance) to r.
     **************************************************************************************************************/

    for (unsigned j = 0; j < my_index_size; j++) {
        unsigned barcode_id = index[j * read_count + read_id].barcode_id;
        barcode b = barcodes[barcode_id];
        index[j * read_count + read_id].distance = exact_distance(b, r, UINT8_MAX);
    }

    /*****************************************************************************************************
     * Sort the index by their final SL-distance. (Insertion-Sort)
     ****************************************************************************************************/

    for (unsigned i = 1; i < my_index_size; i++) {
        // insert my_index[i] into the sorted index [0...i-1]
        index_element el = index[i * read_count + read_id];
        unsigned j = i;
        for (; j > 0 && index[(j - 1) * read_count + read_id].distance > el.distance; j--)
            index[j * read_count + read_id] = index[(j - 1) * read_count + read_id];
        index[j * read_count + read_id] = el;
    }
}


struct k_mer_filtered_calling_gpu_v1 : public barcode_assignment {

    k_mer_filtered_calling_gpu_v1(
        const std::vector<barcode>& barcodes,
        const std::vector<read>& reads,
        const unsigned int k,
        const int distance_measure = SEQUENCE_LEVENSHTEIN_DISTANCE,
        const int rejection_threshold = REJECTION_THRESHOLD)
        : barcode_assignment(run(barcodes, reads, k, distance_measure, rejection_threshold)) {}

private:

    /**
     * Run the calling algorithm with arbitrary values of k and arbitrary distance measure.
     *
     * @param barcodes
     * @param reads
     * @param k
     * @param distance_measure
     * @param rejection_threshold
     * @return
     */
    barcode_assignment run(const std::vector<barcode>& barcodes,
                              const std::vector<read>& reads,
                              const unsigned int k,
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
            return barcode_assignment(reads.size());
        }
    }

    /**
     * Run the calling algorithm with fixed value of k.
     *
     * @tparam k
     * @param barcodes
     * @param reads
     * @param distance_measure
     * @param rejection_threshold
     * @return
     */
    template <unsigned k>
    barcode_assignment run(
        const std::vector<barcode>& barcodes,
        const std::vector<read>& reads,
        const int distance_measure,
        const int rejection_threshold) {

        switch (distance_measure) {
        case SEQUENCE_LEVENSHTEIN_DISTANCE:
            return run<k, sequence_levenshtein_distance_v5>(barcodes, reads, rejection_threshold);
        case LEVENSHTEIN_DISTANCE:
            return run<k, levenshtein_distance_v2>(barcodes, reads, rejection_threshold);
        default:
            std::cerr << "Unknown distance measure" << std::endl;
            exit(1);
        }
    }


    /**
     * Run the calling algorithm with fixed value of k and fixed distance measure.
     *
     * @tparam k K-mer length
     * @tparam exact_distance Distance measure.
     * @param barcodes
     * @param reads
     * @param rejection_threshold
     * @return
     */
    template <unsigned k, uint8_t (*exact_distance)(const sequence&, const sequence&, uint8_t upper_bound)>
    barcode_assignment run(const std::vector<barcode>& barcodes,
                              const std::vector<read>& reads,
                              const int rejection_threshold) {

        static_assert(k >= 2 && k <= 7, "Invalid value of k!");

        barcode_assignment ass(reads.size());

        /************************************************************************************************
         * Allocate gpu memory for the barcodes, reads, and indices.
         ************************************************************************************************/

        barcode* barcodes_dev; // barcode array on gpu memory
        read* reads_dev; // read array on gpu memory
        index_element* index_elements_dev;
        index_element* index_elements_host;

        cudaMalloc(&barcodes_dev, barcodes.size() * sizeof(barcode));
        cudaMalloc(&reads_dev, reads.size() * sizeof(read));
        cudaMalloc(&index_elements_dev, reads.size() * MAX_INDEX_SIZE * sizeof(index_element));

        /************************************************************************************************
         * Copy the barcodes and reads to gpu memory.
         ************************************************************************************************/

        cudaMemcpy(barcodes_dev, &barcodes[0], barcodes.size() * sizeof(barcode), cudaMemcpyHostToDevice);
        cudaMemcpy(reads_dev, &reads[0], reads.size() * sizeof(read), cudaMemcpyHostToDevice);

        /************************************************************************************************
         * Call the calling kernel.
         ************************************************************************************************/

        unsigned threads_per_block = 256;
        unsigned block_count = (reads.size() + threads_per_block - 1) / threads_per_block;
        calling_kernel_kmer_filter_gpu_v1<k, exact_distance>
            <<<block_count, threads_per_block>>>
            (barcodes_dev, reads_dev, index_elements_dev, barcodes.size(), reads.size());

        auto err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << cudaGetErrorString(err) << std::endl;
            return ass;
        }

        /************************************************************************************************
         * Copy the barcode indices to the host memory and select the first of each barcode.
         ************************************************************************************************/

        index_elements_host = new index_element[reads.size() * MAX_INDEX_SIZE];
        cudaMemcpy(index_elements_host, index_elements_dev, reads.size() * MAX_INDEX_SIZE * sizeof(index_element),
                   cudaMemcpyDeviceToHost);

        for (unsigned read_id = 0; read_id < reads.size(); read_id++) {
            if (index_elements_host[read_id].distance <= rejection_threshold)
                ass[read_id] = index_elements_host[read_id].barcode_id;
        }

        /************************************************************************************************
         * Output the index.
         ************************************************************************************************/

        /*for (unsigned read_id = 0; read_id < reads.size(); read_id++) {
            std::cout << "index for read " << read_id << ":" << std::endl;
            for (unsigned i = 0; i < MAX_INDEX_SIZE; i++)
            {
                std::cout << i << " " << index_elements_host[i * reads.size() + read_id].barcode_id << " "
                    << index_elements_host[i * reads.size() + read_id].distance << std::endl;
            }
        }*/

        /************************************************************************************************
         * Free the memory.
         ************************************************************************************************/

        cudaFree(barcodes_dev);
        cudaFree(reads_dev);
        cudaFree(index_elements_dev);
        delete[] index_elements_host;

        return ass;
    }
};


#endif //INC_2OPT_KMER_FILTERED_CALLING_GPU_V1_CUH
