//
// Created by steffen on 24.07.24.
//

#include <iostream>
#include <vector>
#include <chrono>
#include "../precision.h"
#include "k_mer_filtered_calling_host_v1.h"
#include "k_mer_filtered_calling_host_v2.h"
#include "k_mer_filtered_calling_gpu_v4.cuh"
#include "k_mer_filtered_calling_gpu_v5.cuh"
#include "k_mer_filtered_calling_gpu_v6.cuh"
#include "two_step_k_mer_filtered_calling_host_v1.h"
#include "two_step_k_mer_filtered_calling_gpu_v1.cuh"
#include "../recall.h"
#include <omp.h>
#include <thread>
#include "../distance/levenshtein_v1.cuh"
#include "../distance/levenshtein_v2.cuh"
#include "../distance/levenshtein_v3.cuh"
#include "../distance/sequence_levenshtein_v3.cuh"
#include "../distance/sequence_levenshtein_v4.cuh"
#include "../distance/weighted_levenshtein_v1.cuh"
#include "../distance/weighted_sequence_levenshtein_v1.cuh"
#include "../barcode_set_bc_reader.h"
#include "../read_set_fastq_reader.h"
#include "../barcode_assignment_reader.h"
#include "dummy_algorithm_wait_1s.h"

using namespace barcode_calling;

void benchmark_algorithm(const barcode_calling_algorithm& algorithm,
                         const barcode_set& barcodes,
                         const read_set& reads,
                         const barcode_assignment& ground_truth
) {
    // start measuring running time
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration(0);

    // run the algorithm and compute a barcode assignment
    auto ass = algorithm.run(barcodes, reads);

    // stop timing
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;

    std::cout << algorithm.get_name() << "\t"
        << duration.count() / reads.size() << "\t"
        << precision(ass, ground_truth, reads.size()) << "\t"
        << recall(ass, reads.size()) << "\t"
        << std::endl;
}


int main(int argc, char** argv) {

    if (argc != 6) {
        std::cerr <<
            "Usage: ./benchmark_barcode_calling <barcode_file> <read_file> <ground_truth> <distance_measure> <rejection_threshold>"
            << std::endl;
        std::cerr << "  barcode_file: \n"
            "      BC file containing the barcodes (Two lines per barcode)" << std::endl;
        std::cerr << "  read_file: \n"
            "      FASTQ file containing the reads." << std::endl;
        std::cerr << "  ground_truth: \n"
            "      Assignment file in which each read is assigned to its original barcode." << std::endl;
        std::cerr << "  distance_measure: One of the following:\n" <<
            "       levenshtein\n" <<
            "       sequence-levenshtein" <<
            std::endl;
        std::cerr << "  rejection_threshold:\n"
            "       A read is assigned to a barcode only if their associated distance is not\n"
            "       larger than this integer." << std::endl;
        return 0;
    }

    // convert command line arguments
    std::string barcode_file(argv[1]);
    std::string read_file(argv[2]);
    std::string label_file(argv[3]);
    std::string distance_measure_str(argv[4]);
    int rejection_threshold = atoi(argv[5]);

    // read input files from disk
    auto barcodes = barcode_set_bc_reader(barcode_file);
    auto reads = read_set_fastq_reader(read_file);
    auto ground_truth = barcode_assignment_reader(label_file, barcodes, reads);

    // parse distance measure string
    std::unique_ptr<distance_measure> distance;
    if (distance_measure_str == "levenshtein")
        distance = std::make_unique<weighted_levenshtein_v1>(unit_costs());
    else if (distance_measure_str == "sequence-levenshtein")
        distance = std::make_unique<weighted_sequence_levenshtein_v1>(unit_costs());
    else {
        std::cerr << "unsupported distance measure: " << distance_measure_str << std::endl;
        return 5;
    }

    std::cout << "barcode_count               = " << barcodes.size() << std::endl;
    std::cout << "read_count                  = " << reads.size() << std::endl;
    std::cout << "BARCODE_LENGTH              = " << BARCODE_LENGTH << std::endl;
    std::cout << "MAX_INDEX_SIZE              = " << DEFAULT_CANDIDATE_COUNT << std::endl;
    std::cout << "DISTANCE_MEASURE            = " << distance->get_name() << std::endl;
    std::cout << "REJECTION_THRESHOLD         = " << rejection_threshold << std::endl;
    std::cout << "PSEUDO_DISTANCE_WINDOW_SIZE = " << DEFAULT_PSEUDO_DISTANCE_WINDOW_SIZE << std::endl;
    std::cout << "OMP_NUM_THREADS             = " << omp_get_max_threads() << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << "strategy\tms_per_read\tprecision\trecall" << std::endl;


    /****************************************************
     * k-mer-filterted_calling_host_v2
     * k=4,5,6,7
     ***************************************************/

    benchmark_algorithm(k_mer_filtered_calling_host_v2<4>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_host_v2<5>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_host_v2<6>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_host_v2<7>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);

    /****************************************************
     * k-mer-filterted_calling_gpu_v4
     * k=4,5,6,7
     ***************************************************/

    benchmark_algorithm(k_mer_filtered_calling_gpu_v4<4>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_gpu_v4<5>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_gpu_v4<6>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_gpu_v4<7>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);

    /****************************************************
     * k-mer-filterted_calling_gpu_v5
     * k=4,5,6,7
     ***************************************************/

    benchmark_algorithm(k_mer_filtered_calling_gpu_v5<4, 2048>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_gpu_v5<5, 2048>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_gpu_v5<6, 2048>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_gpu_v5<7, 2048>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);

    /****************************************************
     * k-mer-filterted_calling_gpu_v6
     * k=4,5,6,7
     ***************************************************/

    benchmark_algorithm(k_mer_filtered_calling_gpu_v6<4, 512>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_gpu_v6<5, 512>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_gpu_v6<6, 512>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(k_mer_filtered_calling_gpu_v6<7, 512>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);

    /****************************************************
     * two_step_ilterted_calling_host_v1
     * k_large=5,6,7
     * k_small=4
     ***************************************************/

    benchmark_algorithm(two_step_k_mer_filtered_calling_host_v1<7, 4>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(two_step_k_mer_filtered_calling_host_v1<6, 4>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(two_step_k_mer_filtered_calling_host_v1<5, 4>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);

    /****************************************************
     * two_step_ilterted_calling_gpu_v1
     * k_large=5,6,7
     * k_small=4
     ***************************************************/

    benchmark_algorithm(two_step_k_mer_filtered_calling_gpu_v1<7, 4>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(two_step_k_mer_filtered_calling_gpu_v1<6, 4>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);
    benchmark_algorithm(two_step_k_mer_filtered_calling_gpu_v1<5, 4>
                        (*distance, rejection_threshold), barcodes, reads, ground_truth);

    return 0;
}
