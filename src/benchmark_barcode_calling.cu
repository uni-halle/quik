//
// Created by steffen on 24.07.24.
//

#include <iostream>
#include <vector>
#include <chrono>
#include "read_sequence_file.h"
#include "read_label_file.h"
#include "precision.h"
#include "k_mer_filtered_calling_gpu_v1.cuh"
#include "k_mer_filtered_calling_host_v2.h"
#include "k_mer_filtered_calling_gpu_v2.cuh"
#include "k_mer_filtered_calling_gpu_v3.cuh"
#include "k_mer_filtered_calling_gpu_v5.cuh"
#include "k_mer_filtered_calling_gpu_v4.cuh"
#include "recall.h"
#include "two_step_k_mer_filtered_calling_host_v1.h"
#include "two_step_k_mer_filtered_calling_gpu_v1.cuh"
#include <omp.h>
#include <thread>

int main(int argc, char** argv) {

    if (argc != 6) {
        std::cerr << "Usage: ./benchmark <barcode_file> <read_file> <label_file> <distance_measure> <rejection_threshold>"
            << std::endl;
        std::cerr << "  barcode_file: \n"
                     "      Text file with one barcode per line: barcode[0]...barcode[n-1]." << std::endl;
        std::cerr << "  read_file: \n"
                     "      Text file with one read per line: read[0]...read[m-1]." << std::endl;
        std::cerr << "  label_file: \n"
                     "      Text file with m lines: label[0]...label[m-1].\n"
                     "      The label[i] is associated to the read[i] and describes the\n"
                     "      index of the barcode, from which this read originated.\n"
                     "      Thus, read[i] originated from barcode[label[i]]."<< std::endl;
        std::cerr << "  distance_measure: LEVENSHTEIN or SEQUENCE_LEVENSHTEIN" << std::endl;
        return 0;
    }

    // convert command line arguments
    std::string barcode_file(argv[1]);
    std::string read_file(argv[2]);
    std::string label_file(argv[3]);
    std::string distance_measure_str(argv[4]);
    int rejection_threshold = atoi(argv[5]);

    // parse the distance measure
    int distance_measure = -1;
    if (distance_measure_str == "LEVENSHTEIN")
        distance_measure = LEVENSHTEIN_DISTANCE;
    else if (distance_measure_str == "SEQUENCE_LEVENSHTEIN")
        distance_measure = SEQUENCE_LEVENSHTEIN_DISTANCE;
    else {
        std::cerr << "Unsupported distance measure" << std::endl;
        exit(-1);
    }

    // read barcode file
    auto barcodes = read_sequence_file(barcode_file);

    // load read file
    auto reads = read_sequence_file(read_file);

    // load ground truth from the label file
    auto ground_truth = read_label_file(label_file, reads.size());

    // start measuring running time
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration;

    std::vector<barcode_assignment> assignments;

    std::cout << "barcode_count               = " << barcodes.size() << std::endl;
    std::cout << "read_count                  = " << reads.size() << std::endl;
    std::cout << "SEQUENCE_LENGTH             = " << SEQUENCE_LENGTH << std::endl;
    std::cout << "MAX_INDEX_SIZE              = " << MAX_INDEX_SIZE << std::endl;
    std::cout << "DISTANCE_MEASURE            = " << distance_measure_str << std::endl;
    std::cout << "REJECTION_THRESHOLD         = " << rejection_threshold << std::endl;
    std::cout << "PSEUDO_DISTANCE_WINDOW_SIZE = " << PSEUDO_DISTANCE_WINDOW_SIZE << std::endl;
    std::cout << "OMP_NUM_THREADS             = " << omp_get_max_threads() << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << "strategy\tms_per_read\tprecision\trecall" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(k_mer_filtered_calling_gpu_v4(barcodes, reads, 4, distance_measure, rejection_threshold));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "4_mer_filtered_calling_gpu_v4\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(k_mer_filtered_calling_gpu_v4(barcodes, reads, 5, distance_measure, rejection_threshold));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "5_mer_filtered_calling_gpu_v4\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(k_mer_filtered_calling_gpu_v4(barcodes, reads, 6, distance_measure, rejection_threshold));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "6_mer_filtered_calling_gpu_v4\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(k_mer_filtered_calling_gpu_v4(barcodes, reads, 7, distance_measure, rejection_threshold));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "7_mer_filtered_calling_gpu_v4\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(k_mer_filtered_calling_host_v2(barcodes, reads, 4, distance_measure, rejection_threshold));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "4_mer_filtered_calling_host_v2\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(k_mer_filtered_calling_host_v2(barcodes, reads, 5, distance_measure, rejection_threshold));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "5_mer_filtered_calling_host_v2\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(k_mer_filtered_calling_host_v2(barcodes, reads, 6, distance_measure, rejection_threshold));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "6_mer_filtered_calling_host_v2\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(k_mer_filtered_calling_host_v2(barcodes, reads, 7, distance_measure, rejection_threshold));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "7_mer_filtered_calling_host_v2\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    int rejection_threshold_large = 7;
    int rejection_threshold_small = rejection_threshold;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(
        two_step_k_mer_filtered_calling_host_v1(barcodes, reads,
            4, 7, distance_measure, rejection_threshold_small, rejection_threshold_large));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "4_7_mer_filtered_calling_host_v1\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(
        two_step_k_mer_filtered_calling_host_v1(barcodes, reads, 5, 7, distance_measure, rejection_threshold_small, rejection_threshold_large));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "5_7_mer_filtered_calling_host_v1\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(
        two_step_k_mer_filtered_calling_host_v1(barcodes, reads, 6, 7, distance_measure, rejection_threshold_small, rejection_threshold_large));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "6_7_mer_filtered_calling_host_v1\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(
        two_step_k_mer_filtered_calling_gpu_v1(barcodes, reads, 4, 7, distance_measure, rejection_threshold_small, rejection_threshold_large));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "4_7_mer_filtered_calling_gpu_v1\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(
        two_step_k_mer_filtered_calling_gpu_v1(barcodes, reads, 5, 7, distance_measure, rejection_threshold_small, rejection_threshold_large));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "5_7_mer_filtered_calling_gpu_v1\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    assignments.push_back(
        two_step_k_mer_filtered_calling_gpu_v1(barcodes, reads, 6, 7, distance_measure, rejection_threshold_small, rejection_threshold_large));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "6_7_mer_filtered_calling_gpu_v1\t"
        << duration.count() / reads.size() << "\t"
        << precision(assignments.back(), ground_truth) << "\t"
        << recall(assignments.back(), ground_truth) << "\t"
        << std::endl;

    start = std::chrono::high_resolution_clock::now();
    std::this_thread::sleep_for(std::chrono::milliseconds(10000));
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "sleep_10_s" << "\t"
        << duration.count() / reads.size() << "\t"
        << std::endl;

    return 0;
}
