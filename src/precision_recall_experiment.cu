//
// Created by steffen on 02.07.24.
//

#include <iostream>
#include "read_sequence_file.h"
#include "read.h"
#include "exact_calling_gpu_v1.cuh"
#include "exact_calling_gpu_v2.cuh"
#include "exact_calling_gpu_v3.cuh"
#include "exact_calling_gpu_v4.cuh"
#include "exact_calling_gpu_v5.cuh"
#include "k_mer_filtered_calling_gpu_v1.cuh"
#include "k_mer_filtered_calling_gpu_v2.cuh"
#include "k_mer_filtered_calling_gpu_v3.cuh"
#include "k_mer_filtered_calling_host_v1.h"
#include "k_mer_filtered_calling_host_v2.h"
#include "scramble.h"
#include "k_mer_filtered_calling_gpu_v4.cuh"
#include "k_mer_filtered_calling_gpu_v5.cuh"
#include <random>
#include <cassert>
#include <chrono>

inline bool ends_with(std::string const &value, std::string const &ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}


int main(int argc, char **argv) {

    if (argc != 5) {
        std::cerr << "Usage: ./precision_recall_experiment <strategy> <barcode_file> <error_probability> <read_count>"
                  << std::endl;
        std::cerr << " strategy: exact_gpu_v1" << std::endl;
        std::cerr << "           exact_gpu_v2" << std::endl;
        std::cerr << "           exact_gpu_v3" << std::endl;
        std::cerr << "           exact_gpu_v4" << std::endl;
        std::cerr << "           exact_gpu_v5" << std::endl;
        std::cerr << "           k_mer_filter_gpu_v1 (replace k by an integer between 2 and 8)" << std::endl;
        std::cerr << "           k_mer_filter_gpu_v2 (replace k by an integer between 2 and 8)" << std::endl;
        std::cerr << "           k_mer_filter_gpu_v3 (replace k by an integer between 2 and 8)" << std::endl;
        std::cerr << "           k_mer_filter_gpu_v4 (replace k by an integer between 2 and 8)" << std::endl;
        std::cerr << "           k_mer_filter_gpu_v5 (replace k by an integer between 2 and 8)" << std::endl;
        std::cerr << "           k_mer_filter_host_v1 (replace k by an integer between 2 and 8)" << std::endl;
        std::cerr << "           k_mer_filter_host_v2 (replace k by an integer between 2 and 8)" << std::endl;
        return 0;
    }

    // convert command line arguments
    std::string strategy(argv[1]);
    std::string barcode_file(argv[2]);
    double p = atof(argv[3]);
    int read_count = atoi(argv[4]);

    // read barcode file
    auto barcodes = read_sequence_file(barcode_file);

    // simulate sequencing errors
    std::vector<read> reads;
    std::vector<unsigned> labels;
    while (reads.size() < read_count) {
        for (unsigned barcode_id = 0; barcode_id < barcodes.size() && reads.size() < read_count; barcode_id++) {
            labels.push_back(barcode_id);
            reads.emplace_back(scramble(barcodes[barcode_id], p / 3, p / 3, p / 3));
        }
    }

    // start measuring running time
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<unsigned> assignment;

    // apply barcode calling scheme
    if (strategy == "exact_gpu_v1") {
        assignment = exact_calling_gpu_v1(barcodes, reads);
    } else if (strategy == "exact_gpu_v2") {
        assignment = exact_calling_gpu_v2(barcodes, reads);
    } else if (strategy == "exact_gpu_v3") {
        assignment = exact_calling_gpu_v3(barcodes, reads);
    } else if (strategy == "exact_gpu_v4") {
        assignment = exact_calling_gpu_v4(barcodes, reads);
    } else if (strategy == "exact_gpu_v5") {
        assignment = exact_calling_gpu_v5(barcodes, reads);
    } else if (ends_with(strategy, "_mer_filter_gpu_v1")) {
        int k = std::stoi(strategy.substr(0, 1));
        assignment = k_mer_filtered_calling_gpu_v1(barcodes, reads, k);
    } else if (ends_with(strategy, "_mer_filter_gpu_v2")) {
        int k = std::stoi(strategy.substr(0, 1));
        assignment = k_mer_filtered_calling_gpu_v2(barcodes, reads, k);
    } else if (ends_with(strategy, "_mer_filter_gpu_v3")) {
        int k = std::stoi(strategy.substr(0, 1));
        assignment = k_mer_filtered_calling_gpu_v3(barcodes, reads, k);
    } else if (ends_with(strategy, "_mer_filter_gpu_v4")) {
        int k = std::stoi(strategy.substr(0, 1));
        assignment = k_mer_filtered_calling_gpu_v4(barcodes, reads, k);
    } else if (ends_with(strategy, "_mer_filter_gpu_v5")) {
        int k = std::stoi(strategy.substr(0, 1));
        assignment = k_mer_filtered_calling_gpu_v5(barcodes, reads, k);
    } else if (ends_with(strategy, "_mer_filter_host_v1")) {
        int k = std::stoi(strategy.substr(0, 1));
        assignment = k_mer_filtered_calling_host_v1(barcodes, reads, k);
    } else if (ends_with(strategy, "_mer_filter_host_v2")) {
        int k = std::stoi(strategy.substr(0, 1));
        assignment = k_mer_filtered_calling_host_v2(barcodes, reads, k);
    } else {
        std::cerr << "Error! Unknown strategy: " << strategy << std::endl;
        return -2;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;

    // calculate precision and recall
    unsigned accepted_reads = 0;
    unsigned correctly_identified_reads = 0;
    for (unsigned read_id = 0; read_id < reads.size(); read_id++) {
        if (assignment[read_id] != (unsigned) -1) {
            accepted_reads++;
            if (assignment[read_id] == labels[read_id])
                correctly_identified_reads++;
        }
    }

    double precision = (double) correctly_identified_reads / accepted_reads;
    double recall = (double) accepted_reads / reads.size();

    std::cout << "strategy       = " << strategy << std::endl;
    std::cout << "barcode_count  = " << barcodes.size() << std::endl;
    std::cout << "read_count     = " << read_count << std::endl;
    std::cout << "MAX_INDEX_SIZE = " << MAX_INDEX_SIZE << std::endl;
    std::cout << "accepted_reads = " << accepted_reads << std::endl;
    std::cout << "correct_reads  = " << correctly_identified_reads << std::endl;
    std::cout << "precision      = " << precision << std::endl;
    std::cout << "recall         = " << recall << std::endl;
    std::cout << "ms per read    = " << duration.count() / (double) reads.size() << std::endl;

    /*std::cout << "not correctly identified reads=" << std::endl;
    for (unsigned read_id = 0; read_id < reads.size(); read_id++) {
        if (assignment[read_id] != (unsigned) -1) {
            if (assignment[read_id] != labels[read_id])
                std::cout << read_id << std::endl;
        }
    }*/

    return 0;
}