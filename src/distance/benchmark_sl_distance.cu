//
// Created by steffen on 24.07.24.
//

#include <iostream>
#include <vector>
#include <chrono>
#include "read_sequence_file.h"
#include "scramble.h"
#include "sequence_levenshtein_v1.cuh"
#include "sequence_levenshtein_distance_v2.cuh"
#include "sequence_levenshtein_distance_v3.cuh"
#include "sequence_levenshtein_distance_v4.cuh"
#include "sequence_levenshtein_distance_v5.cuh"
#include <omp.h>
#include <thread>

template <
    uint8_t (*distance)(const sequence&, const sequence&, uint8_t upper_bound) // distance measure
>
__global__ void benchmark_sl_distance_kernel(const barcode* barcodes,
                                             const read* reads,
                                             uint8_t* min_distances,
                                             const size_t barcode_cnt,
                                             const size_t read_cnt) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    // for each read
    for (size_t read_id = id; read_id < read_cnt; read_id += stride) {
        uint8_t min_sl_dist = UINT8_MAX;
        read r = reads[read_id];
        // for each barcode
        for (size_t barcode_id = 0; barcode_id < barcode_cnt; barcode_id++) {
            barcode b = barcodes[barcode_id];
            // calculate the sl-distance between b and r
            uint8_t dist = distance(b, r, min_sl_dist);
            if (dist < min_sl_dist)
                min_sl_dist = dist;
        }

        min_distances[read_id] = min_sl_dist;
    }
}

template <
    uint8_t (*distance)(const sequence&, const sequence&, uint8_t upper_bound) // distance measure
>
double benchmark_sl_distance(const std::vector<sequence>& barcodes, const std::vector<sequence>& reads) {

    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration;

    /************************************************************************************************
   * Allocate gpu memory for barcodes and reads.
   ************************************************************************************************/

    barcode* barcodes_dev; // barcode array on gpu memory
    read* reads_dev; // read array on gpu memory
    uint8_t* min_distances_dev, *min_distances_host; // minimum distance of each read

    cudaMalloc(&barcodes_dev, barcodes.size() * sizeof(barcode));
    cudaMalloc(&reads_dev, reads.size() * sizeof(read));
    cudaMalloc(&min_distances_dev, reads.size() * sizeof(read));
    min_distances_host = (uint8_t*)malloc(reads.size() * sizeof(uint8_t));

    /************************************************************************************************
     * Copy the barcodes and reads to gpu memory.
     ************************************************************************************************/

    cudaMemcpy(barcodes_dev, &barcodes[0], barcodes.size() * sizeof(barcode), cudaMemcpyHostToDevice);
    cudaMemcpy(reads_dev, &reads[0], reads.size() * sizeof(read), cudaMemcpyHostToDevice);

    /************************************************************************************************
     * Call the benchmark kernel.
     ************************************************************************************************/

    start = std::chrono::high_resolution_clock::now();

    const int block_cnt = 512;
    const int threads_per_block = 256;
    benchmark_sl_distance_kernel
    <distance>
    <<<block_cnt, threads_per_block>>>
    (barcodes_dev, reads_dev, min_distances_dev, barcodes.size(), reads.size());

    cudaDeviceSynchronize();

    end = std::chrono::high_resolution_clock::now();

    /************************************************************************************************
     * Copy back the minimum distances.
     ************************************************************************************************/

    cudaMemcpy(min_distances_host, min_distances_dev, reads.size() * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    /************************************************************************************************
     * Free the memory.
     ************************************************************************************************/

    cudaFree(min_distances_dev);
    cudaFree(barcodes_dev);
    cudaFree(reads_dev);
    free(min_distances_host);

    duration = end - start;

    return duration.count();
}


int main(int argc, char** argv) {

    if (argc != 3) {
        std::cerr << "Usage: ./benchmark_sl_distance <barcode_file> <read_file>"
            << std::endl;
        return 0;
    }

    // convert command line arguments
    std::string barcode_file(argv[1]);
    std::string read_file(argv[2]);

    // read barcode file
    auto barcodes = read_sequence_file(barcode_file);

    // load read file
    auto reads = read_sequence_file(read_file);

    // start measuring running time


    std::cout << "barcode_count               = " << barcodes.size() << std::endl;
    std::cout << "read_count                  = " << reads.size() << std::endl;
    std::cout << "SEQUENCE_LENGTH             = " << SEQUENCE_LENGTH << std::endl;
    std::cout << "OMP_NUM_THREADS             = " << omp_get_max_threads() << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << "strategy\tms_per_read" << std::endl;

    std::cout << "sequence_levenshtein_v1" << "\t"
        << benchmark_sl_distance<sequence_levenshtein_v1>(barcodes, reads) / reads.size() << "\t"
        << std::endl;

    std::cout << "sequence_levenshtein_distance_v2" << "\t"
        << benchmark_sl_distance<sequence_levenshtein_distance_v2>(barcodes, reads) / reads.size() << "\t"
        << std::endl;

    std::cout << "sequence_levenshtein_distance_v3" << "\t"
        << benchmark_sl_distance<sequence_levenshtein_distance_v3>(barcodes, reads) / reads.size() << "\t"
        << std::endl;

    std::cout << "sequence_levenshtein_distance_v4" << "\t"
        << benchmark_sl_distance<sequence_levenshtein_distance_v4>(barcodes, reads) / reads.size() << "\t"
        << std::endl;

    std::cout << "sequence_levenshtein_distance_v5" << "\t"
        << benchmark_sl_distance<sequence_levenshtein_distance_v5>(barcodes, reads) / reads.size() << "\t"
        << std::endl;

    return 0;
}
