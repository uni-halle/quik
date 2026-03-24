//
// Created by steffen on 24.07.24.
//

#include <iostream>
#include <vector>
#include <chrono>
#include "../precision.h"
#include "k_mer_filter_host_v2.h"
#include "k_mer_filter_single_gpu_v09.cuh"
#include "../recall.h"
#include <omp.h>
#include <thread>
#include "constants.h"
#include "k_mer_filter_multi_gpu_default.cuh"
#include "two_step_k_mer_filter_host_v1.h"
#include "two_step_k_mer_filter_single_gpu_v1.cuh"
#include "two_step_k_mer_filter_single_gpu_v2.cuh"
#include "two_step_k_mer_filter_multi_gpu_v1.cuh"
#include "two_step_k_mer_filter_multi_gpu_v2.cuh"
#include "../distance/weighted_levenshtein_v1.cuh"
#include "../distance/weighted_sequence_levenshtein_v1.cuh"
#include "../barcode_set_bc_reader.h"
#include "../read_file_fastq_reader.h"
#include "../barcode_assignment_reader.h"

using namespace barcode_calling;

template <typename algorithm>
void benchmark_algorithm(const barcode_set& barcodes,
                         const read_set& reads,
                         const int rejection_threshold,
                         const barcode_assignment& ground_truth) {

    // create the algorithm object
    algorithm alg(barcodes, reads, rejection_threshold);

    int trials = 3;
    double seconds_total = 0;

    for (int i=0; i<trials; i++) {

        // start measuring running time
        auto start = std::chrono::high_resolution_clock::now();

        // run the algorithm and compute a barcode assignment
        alg.run();

        // stop timing
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = end - start;
        auto seconds = std::chrono::duration<double>(elapsed).count();
        seconds_total += seconds;
        //std::cout << i << ": " << seconds << " seconds" << std::endl;
    }

    double seconds = seconds_total / trials;

    double reads_per_sec = reads.size() / seconds;
    double ms_per_read = seconds * 1000 / reads.size();

    std::cout << algorithm::name() << "\t"
        << seconds << "\t"
        << ms_per_read << "\t"
        << reads_per_sec << "\t"
        << precision(reads, alg, ground_truth, rejection_threshold) << "\t"
        << recall(reads, alg, rejection_threshold) << "\t"
        << std::endl;
}

template <typename distance_function>
int run(const barcode_set& barcodes,
        const read_set& reads,
        const int32_t rejection_threshold,
        const barcode_assignment& ground_truth
) {

    // read GPU data
    int deviceCount = 0;
    CUDA_CHECK(cudaGetDeviceCount(&deviceCount));

    std::cout << "barcode_count               = " << barcodes.size() << std::endl;
    std::cout << "read_count                  = " << reads.size() << std::endl;
    std::cout << "BARCODE_LENGTH              = " << BARCODE_LENGTH << std::endl;
    std::cout << "MAX_INDEX_SIZE              = " << DEFAULT_CANDIDATE_COUNT << std::endl;
    std::cout << "DISTANCE_MEASURE            = " << distance_function::name() << std::endl;
    std::cout << "REJECTION_THRESHOLD         = " << rejection_threshold << std::endl;
    std::cout << "PSEUDO_DISTANCE_WINDOW_SIZE = " << DEFAULT_PSEUDO_DISTANCE_WINDOW_SIZE << std::endl;
    std::cout << "OMP_NUM_THREADS             = " << omp_get_max_threads() << std::endl;
    std::cout << "GPU count                   = " << deviceCount << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << "strategy\tsec\tms_per_read\treads_per_sec\tprecision\trecall" << std::endl;

    /****************************************************
     * 4-mer-filter
     ***************************************************/

    benchmark_algorithm<k_mer_filter_host_v2<4, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<k_mer_filter_single_gpu_v09<4, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<k_mer_filter_multi_gpu_default<4, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    /****************************************************
     * 5-mer-filter
     ***************************************************/

    benchmark_algorithm<k_mer_filter_host_v2<5, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<k_mer_filter_single_gpu_v09<5, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<k_mer_filter_multi_gpu_default<5, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    /****************************************************
     * 6-mer-filter
     ***************************************************/

    benchmark_algorithm<k_mer_filter_host_v2<6, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<k_mer_filter_single_gpu_v09<6, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<k_mer_filter_multi_gpu_default<6, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    /****************************************************
     * 7-mer-filter
     ***************************************************/

    benchmark_algorithm<k_mer_filter_host_v2<7, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<k_mer_filter_single_gpu_v09<7, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<k_mer_filter_multi_gpu_default<7, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    /****************************************************
    * <7,4>-mer-filter
    ***************************************************/

    benchmark_algorithm<two_step_k_mer_filter_host_v1<7, 4, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<two_step_k_mer_filter_single_gpu_v1<7, 4, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<two_step_k_mer_filter_single_gpu_v2<7, 4, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<two_step_k_mer_filter_multi_gpu_v1<7, 4, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<two_step_k_mer_filter_multi_gpu_v2<7, 4, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);


    /****************************************************
     * <7,5>-mer-filter
     ***************************************************/

    benchmark_algorithm<two_step_k_mer_filter_host_v1<7, 5, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<two_step_k_mer_filter_single_gpu_v1<7, 5, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<two_step_k_mer_filter_single_gpu_v2<7, 5, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<two_step_k_mer_filter_multi_gpu_v1<7, 5, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    benchmark_algorithm<two_step_k_mer_filter_multi_gpu_v2<7, 5, distance_function>>
        (barcodes, reads, rejection_threshold, ground_truth);

    return 0;
}

int main(int argc, char** argv) {

    if (argc != 6) {
        std::cerr <<
            "Usage: ./benchmark_barcode_calling <barcode_file> <read_file> <label_file> <distance_measure> <rejection_threshold>"
            << std::endl;
        std::cerr << "  barcode_file: \n"
            "      Text file with one barcode per line: barcode[0]...barcode[n-1]." << std::endl;
        std::cerr << "  read_file: \n"
            "      Text file with one read per line: read[0]...read[m-1]." << std::endl;
        std::cerr << "  label_file: \n"
            "      Text file with m lines: label[0]...label[m-1].\n"
            "      The label[i] is associated to the read[i] and describes the\n"
            "      index of the barcode, from which this read originated.\n"
            "      Thus, read[i] originated from barcode[label[i]]." << std::endl;
        std::cerr << "  distance_measure: One of the following:\n" <<
            "       levenshtein\n" <<
            "       sequence_levenshtein\n" <<
            std::endl;
        std::cerr << "  rejection_threshold: "
            "       A read is assigned to the closest barcodes only if its distance not "
            "       larger than the rejection threshold." << std::endl;
        return 0;
    }

    // convert command line arguments
    std::string barcode_file_str(argv[1]);
    std::string read_file_str(argv[2]);
    std::string label_file_str(argv[3]);
    std::string distance_measure_str(argv[4]);
    int rejection_threshold = atoi(argv[5]);

    // read input files from disk
    auto barcodes = barcode_set_bc_reader(barcode_file_str);
    auto read_file = read_file_fastq_reader(read_file_str).next();
    auto reads = read_set(read_file);
    auto ground_truth = barcode_assignment_reader(label_file_str, barcodes, read_file);

    // parse distance measure string
    if (distance_measure_str == "levenshtein")
        return run<weighted_levenshtein_v1<unit_costs>>(
            barcodes, reads, rejection_threshold, ground_truth);

    if (distance_measure_str == "sequence_levenshtein")
        return run<weighted_sequence_levenshtein_v1<unit_costs>>(
            barcodes, reads, rejection_threshold, ground_truth);

    std::cerr << "unsupported distance measure: " << distance_measure_str << std::endl;
    return 5;
}
