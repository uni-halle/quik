//
// Created by steffen on 10.07.25.
//

#include <iostream>
#include <vector>
#include <chrono>
#include "calling_algorithms/k_mer_filtered_calling_host_v2.h"
#include <omp.h>
#include <thread>
#include <unordered_map>

#include "barcode_assignment_writer.h"
#include "barcode_set_bc_reader.h"
#include "read_set_fastq_reader.h"
#include "distance/levenshtein_v3.cuh"
#include "distance/weighted_levenshtein_v1.cuh"
#include "distance/sequence_levenshtein_v3.cuh"
#include "extended_barcode_assignment_writer.h"
#include "number_of_assigned_reads.h"
#include "calling_algorithms/k_mer_filtered_calling_gpu_v5.cuh"
#include "calling_algorithms/two_step_k_mer_filtered_calling_gpu_v1.cuh"
#include "calling_algorithms/two_step_k_mer_filtered_calling_host_v1.h"
#include "distance/weighted_sequence_levenshtein_v1.cuh"

using namespace barcode_calling;


void print_help(std::ostream& os = std::cout) {
    os <<
        R"(quik: Fast barcode matching for sequencing reads

Usage:
  quik --barcodes <FILE> --reads <FILE> [OPTIONS]

Required arguments:
  -b, --barcodes <FILE>     Barcode file in our custom BC format.

                            Such files have two line-separated fields per barcode.
                              - Field 1 begins with a '@' character and is followed
                                by a barcode identifier (similar to FASTQ).
                              - Field 2 is the raw sequence letters.

  -r, --reads <FILE>        Read file in FASTQ format.

Optional arguments:
  -d, --distance <STRING>   Distance measure specification

                            Possible values:

                              levenshtein
                              sequence-levenshtein

                            (default: sequence-levenshtein)

  -t, --threshold-distance  Maximum allowed distance between read and barcode.
                            A read is only assigned to a barcode if its
                            distance is smaller or equal than this integer.

                            (default: 2.147.483.647 (INT32_MAX))

  -m, --method <STRING>     Barcode matching method
                            Possible values:

                              4-mer-filter
                                  Highest accuracy, slowest filter method

                              5-mer-filter
                                  Decent accuracy, faster than 4-mer-filter

                              6-mer-filter
                                  Low accuracy, faster than 5-mer-filter

                              7-mer-filter
                                  Lowest accuracy, fastest method

                              7-4-mer-filter
                                  First runs the fast 7-mer-filter to assign the
                                  majority of the reads, then run the 4-mer-filter
                                  to assign the remaining reads.

                            (default: 4-mer-filter)

  -g, --gpu                 Apply all calculations on the first GPU visible to quik.
                            This is usually much faster than non-gpu mode.

  -h, --help                Show this help message and exit
  -v, --verbose             Print extra information to the standard error stream

Output:
  quik writes one line per assigned read to the standard output.
  Fields are tab-separated and have the following meaning:

    read           Sequence identifier of each read (at most once)
    barcode        Sequence identifier of closest barcode
    distance       Distance between read and barcode

  Reads with distance > threshold-distance do not occur in the output.

Examples:
  # Basic usage on GPU
  quik --barcodes barcodes.bc --reads reads.fq --gpu

  # Use 5-mer filtering and classical Levenshtein distance
  quik -b barcodes.bc -r reads.fq -m 5-mer-filter -d levenshtein
)";
}

std::unordered_map<std::string, std::string>
parse_arguments(int argc, char* argv[]) {
    std::unordered_map<std::string, std::string> args;

    for (int i = 1; i < argc; ++i) {
        std::string token = argv[i];

        // Long option: --key or --key=value
        if (token.rfind("--", 0) == 0) {
            auto eq = token.find('=');

            if (eq != std::string::npos) {
                // --key=value
                std::string key = token.substr(2, eq - 2);
                std::string value = token.substr(eq + 1);
                args[key] = value;
            } else {
                // --key [value]  or flag
                std::string key = token.substr(2);

                if (i + 1 < argc && argv[i + 1][0] != '-') {
                    args[key] = argv[++i];
                } else {
                    args[key] = ""; // flag
                }
            }
        }
        // Short option: -k value
        else if (token.rfind('-', 0) == 0) {
            std::string key = token.substr(1);

            if (key.empty())
                throw std::runtime_error("Invalid option '-'");

            if (i + 1 < argc && argv[i + 1][0] != '-') {
                args[key] = argv[++i];
            } else {
                args[key] = ""; // flag
            }
        }
        // Positional argument (optional)
        else {
            args["_positional"].append(token + " ");
        }
    }

    return args;
}

std::string get_arg(const std::unordered_map<std::string, std::string>& args,
                    const std::string& long_opt,
                    const std::string& short_opt,
                    const std::string& default_val = "") {
    if (args.count(long_opt))
        return args.at(long_opt);
    if (args.count(short_opt))
        return args.at(short_opt);
    return default_val;
}


int main(int argc, char** argv) {

    // parse arguments
    auto args = parse_arguments(argc, argv);

    if (args.contains("help") || args.contains("h")) {
        print_help();
        return -1;
    }

    bool gpu = args.count("gpu") || args.count("g");
    bool verbose = args.count("verbose") || args.count("v");

    std::string barcode_file = args.contains("barcodes") ? args["barcodes"] : args["b"];
    std::string read_file = args.contains("reads") ? args["reads"] : args["r"];

    if (barcode_file.empty() || read_file.empty()) {
        std::cerr << "Error! Missing required arguments --barcodes and --reads" << std::endl;
        return 4;
    }

    std::string method_str = get_arg(args, "method", "m", "4-mer-filter");
    std::string distance_str = get_arg(args, "distance", "d", "sequence-levenshtein");
    std::string threshold_str = get_arg(args, "threshold_distance", "t", "");

    // load barcode and reads from disk
    auto barcodes = barcode_set_bc_reader(barcode_file);
    auto reads = read_set_fastq_reader(read_file);

    if (barcodes.empty()) {
        std::cerr << "Error! Empty barcode file!" << std::endl;
        return 1;
    }

    if (reads.empty()) {
        std::cerr << "Error! Empty read file!" << std::endl;
        return 2;
    }

    // read threshold
    int32_t rejection_threshold = INT32_MAX;
    if (!threshold_str.empty())
        rejection_threshold = std::stoi(threshold_str);


    /****************************************************************************
     * Select the distance function. We use weighted variants in all cases.
     ***************************************************************************/

    std::unique_ptr<distance_measure> dist;
    if (distance_str == "levenshtein")
        dist = std::make_unique<weighted_levenshtein_v1>(unit_costs());
    else if (distance_str == "sequence-levenshtein")
        dist = std::make_unique<weighted_sequence_levenshtein_v1>(unit_costs());
    else {
        std::cerr << "Error! Unknown distance!" << std::endl;
        return 3;
    }

    /****************************************************************************
    * Select the calling algorithm.
    ****************************************************************************/

    std::unique_ptr<barcode_calling_algorithm> alg;

    std::string gpu_name = "none";

    if (gpu) {

        /******************************************************************
         * Read basic gpu data
         *****************************************************************/

        // read GPU data
        int deviceCount = 0;
        cudaError_t err = cudaGetDeviceCount(&deviceCount);
        if (err != cudaSuccess) {
            std::cerr << "Error while reading CUDA device count: "
                      << cudaGetErrorString(err) << std::endl;
            return 1;
        }

        if (deviceCount == 0) {
            std::cout << "Error! No CUDA-capable GPUs found!" << std::endl;
            return -4;
        }

        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0); // 0 = first GPU
        gpu_name = prop.name;

        if (method_str == "4-mer-filter")
            alg = std::make_unique<k_mer_filtered_calling_gpu_v5<4>>(*dist);
        else if (method_str == "5-mer-filter")
            alg = std::make_unique<k_mer_filtered_calling_gpu_v5<5>>(*dist);
        else if (method_str == "6-mer-filter")
            alg = std::make_unique<k_mer_filtered_calling_gpu_v5<6>>(*dist);
        else if (method_str == "7-mer-filter")
            alg = std::make_unique<k_mer_filtered_calling_gpu_v5<7>>(*dist);
        else if (method_str == "7-4-mer-filter")
            alg = std::make_unique<two_step_k_mer_filtered_calling_gpu_v1<7,4>>(*dist, rejection_threshold);
        else {
            throw std::runtime_error("Unknown barcode calling method.");
        }
    }
    else {
        if (method_str == "4-mer-filter")
            alg = std::make_unique<k_mer_filtered_calling_host_v2<4>>(*dist);
        else if (method_str == "5-mer-filter")
            alg = std::make_unique<k_mer_filtered_calling_host_v2<5>>(*dist);
        else if (method_str == "6-mer-filter")
            alg = std::make_unique<k_mer_filtered_calling_host_v2<6>>(*dist);
        else if (method_str == "7-mer-filter")
            alg = std::make_unique<k_mer_filtered_calling_host_v2<7>>(*dist);
        else if (method_str == "7-4-mer-filter")
            alg = std::make_unique<two_step_k_mer_filtered_calling_host_v1<7,4>>(*dist, rejection_threshold);
        else {
            throw std::runtime_error("Unknown barcode calling method.");
        }
    }

    // start measuring running time
    auto start = std::chrono::high_resolution_clock::now();

    // run the algorithm to compute a barcode assignment
    auto ass = alg->run(barcodes, reads);

    // stop measuring running time
    auto end = std::chrono::high_resolution_clock::now();

    // output to the standard output stream
    barcode_assignment_writer(barcodes, reads, ass, rejection_threshold).write(std::cout);

    double running_time_sec = std::chrono::duration<double>(end - start).count();
    if (verbose) {

        // print debug information
        std::cerr << "barcode_file             = " << barcode_file << std::endl;
        std::cerr << "read_file                = " << read_file << std::endl;
        std::cerr << "barcode_count            = " << barcodes.size() << std::endl;
        std::cerr << "read_count               = " << reads.size() << std::endl;
        std::cerr << "BARCODE_LENGTH           = " << BARCODE_LENGTH << std::endl;
        std::cerr << "method                   = " << alg->get_name() << std::endl;
        std::cerr << "distance                 = " << dist->get_name() << std::endl;
        std::cerr << "rejection_threshold      = " << rejection_threshold << std::endl;
        std::cerr << "OMP_NUM_THREADS          = " << omp_get_max_threads() << std::endl;
        std::cerr << "GPU                      = " << gpu_name << std::endl;

        unsigned num_assigned = number_of_assigned_reads(reads, ass, rejection_threshold);
        float percent_assigned = 100.0 * num_assigned / reads.size();
        std::cerr << "total running time (sec) = " << running_time_sec << std::endl;
        std::cerr << "number of assigned reads = " << num_assigned << " ("
        << percent_assigned << " %)" << std::endl;
    }

    return 0;
}
