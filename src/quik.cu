//
// Created by steffen on 10.07.25.
//

#include <iostream>
#include <vector>
#include <chrono>
#include "calling_algorithms/k_mer_filter_host_v2.h"
#include <omp.h>
#include <thread>
#include <unordered_map>
#include "explicit_costs.h"
#include "barcode_set_bc_reader.h"
#include "read_file_fastq_reader.h"
#include "distance/weighted_levenshtein_v1.cuh"
#include "extended_barcode_assignment_writer.h"
#include "constants.h"
#include "calling_algorithms/k_mer_filter_host_default.cuh"
#include "calling_algorithms/k_mer_filter_multi_gpu_default.cuh"
#include "calling_algorithms/two_step_k_mer_filter_host_default.cuh"
#include "calling_algorithms/two_step_k_mer_filter_multi_gpu_default.cuh"
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
                                  High accuracy, decent running time

                              5-mer-filter
                                  Decent accuracy, faster than the 4-mer filter

                              6-mer-filter
                                  Low accuracy, faster than the 5-mer filter

                              7-mer-filter
                                  Lowest accuracy, faster than the 6-mer filter

                              7-4-mer-filter
                                  First runs the fast 7-mer-filter to assign the
                                  majority of the reads, then run the 4-mer-filter
                                  to assign the remaining reads.

                            (default: 4-mer-filter)

  -g, --gpu                 Apply the calculations on all CUDA capable GPUs that
                            are available at this system.
                            This is usually much faster than non-gpu mode.

  -h, --help                Show this help message and exit
  -v, --verbose             Print extra information to the standard error stream

Output:

  quik writes one line per assigned read to the standard output.
  Fields are tab-separated and have the following meaning:

    read           Sequence identifier of each read (at most once)
    barcode1       Sequence identifier of closest barcode
    barcode2       Sequence identifier of second to closest barcode
    distance1      Distance between read and barcode1
    distance2      Distance between read and barcode2

  Reads with distance1 > threshold-distance do not occur in the output.

Examples:
  # Basic usage on GPU
  quik --barcodes barcodes.bc --reads reads.fq --gpu

  # Use 5-mer filtering and classical Levenshtein distance
  quik -b barcodes.bc -r reads.fq -m 5-mer-filter -d levenshtein
)";
}

/*************************************************************
 * Command line parameters
 ************************************************************/

std::string barcode_file_str;
std::string read_file_str;
std::string cost_file_str;
std::string threshold_str;
std::string distance_str;
std::string method_str;
bool verbose, gpu;


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

template <typename calling_algorithm, typename distance_function>
int run2(const barcode_set& barcodes,
         read_file_fastq_reader& read_file_reader,
         const int32_t rejection_threshold,
         const distance_function dist) {

    // we count the total number of assigned reads
    size_t assigned_reads = 0;

    // start measuring running time
    auto start = std::chrono::high_resolution_clock::now();

    /*****************************************************************************************
     * While there are still reads in the file, load and process the next chunk of reads
     *****************************************************************************************/

    for (read_file read_f = read_file_reader.next(MAX_READ_COUNT); !read_f.empty();
        read_f = read_file_reader.next(MAX_READ_COUNT)) {

        //std::cout << "process chunk of " << read_f.size() << " reads" << std::endl;

        // convert the read file to the read set
        read_set reads(read_f);

        // create the algorithm object
        calling_algorithm alg(barcodes, reads, rejection_threshold, dist);

        // run the algorithm to compute a barcode assignment
        const extended_barcode_assignment& ass = alg.run();

        // output to the standard output stream
        extended_barcode_assignment_writer(barcodes, read_f, ass, rejection_threshold).write(std::cout);

        // update the assignment count
        for (size_t read_id= 0; read_id < reads.size(); read_id++) {
            if (ass.get_1st_distances()[read_id] <= rejection_threshold)
                assigned_reads++;
        }
    }

    // stop measuring running time
    auto end = std::chrono::high_resolution_clock::now();
    double running_time_sec = std::chrono::duration<double>(end - start).count();

    if (verbose) {

        int device_count = 0;
        CUDA_CHECK(cudaGetDeviceCount(&device_count));

        size_t read_count = read_file_reader.get_processed_read_count();

        // print debug information
        std::cerr << "barcode_file             = " << barcode_file_str << std::endl;
        std::cerr << "read_file                = " << read_file_str << std::endl;
        std::cerr << "barcode_count            = " << barcodes.size() << std::endl;
        std::cerr << "read_count               = " << read_count << std::endl;
        std::cerr << "BARCODE_LENGTH           = " << BARCODE_LENGTH << std::endl;
        std::cerr << "costs                    = " << dist.get_costs().to_string() << std::endl;
        std::cerr << "method                   = " << calling_algorithm::name() << std::endl;
        std::cerr << "distance                 = " << distance_function::name() << std::endl;
        std::cerr << "rejection_threshold      = " << rejection_threshold << std::endl;
        std::cerr << "OMP_NUM_THREADS          = " << omp_get_max_threads() << std::endl;
        std::cerr << "GPU count                = " << device_count << std::endl;

        double percent_assigned = 100.0 * assigned_reads / read_count;
        std::cerr << "total running time (sec) = " << running_time_sec << std::endl;
        std::cerr << "number of assigned reads = " << assigned_reads << " ("
            << percent_assigned << " %)" << std::endl;
    }

    return 0;

}


template <typename distance_function>
int run1(const barcode_set& barcodes,
         read_file_fastq_reader& read_file_reader,
         const int32_t rejection_threshold,
         const distance_function dist) {

    /****************************************************************************
    * Select the calling algorithm.
    ****************************************************************************/

    if (gpu) {

        /******************************************************************
         * Read basic gpu data
         *****************************************************************/

        if (method_str == "4-mer-filter")
            return run2<k_mer_filter_multi_gpu_default<4, distance_function>, distance_function>(
                barcodes, read_file_reader, rejection_threshold, dist);

        if (method_str == "5-mer-filter")
            return run2<k_mer_filter_multi_gpu_default<5, distance_function>, distance_function>(
                barcodes, read_file_reader, rejection_threshold, dist);

        if (method_str == "6-mer-filter")
            return run2<k_mer_filter_multi_gpu_default<6, distance_function>, distance_function>(
                barcodes, read_file_reader, rejection_threshold, dist);

        if (method_str == "7-mer-filter")
            return run2<k_mer_filter_multi_gpu_default<7, distance_function>, distance_function>(
                barcodes, read_file_reader, rejection_threshold, dist);

        if (method_str == "7-4-mer-filter")
            return run2<two_step_k_mer_filter_multi_gpu_default<7, 4, distance_function>, distance_function>(
                barcodes, read_file_reader, rejection_threshold, dist);


        std::cerr << "Error! Unknown barcode calling method" << method_str << std::endl;
        return 7;

    }

    // if no gpu flag was given

    if (method_str == "4-mer-filter")
        return run2<k_mer_filter_host_default<4, distance_function>, distance_function>(
            barcodes, read_file_reader, rejection_threshold, dist);

    if (method_str == "5-mer-filter")
        return run2<k_mer_filter_host_default<5, distance_function>, distance_function>(
            barcodes, read_file_reader, rejection_threshold, dist);

    if (method_str == "6-mer-filter")
        return run2<k_mer_filter_host_default<6, distance_function>, distance_function>(
            barcodes, read_file_reader, rejection_threshold, dist);

    if (method_str == "7-mer-filter")
        return run2<k_mer_filter_host_default<7, distance_function>, distance_function>(
            barcodes, read_file_reader, rejection_threshold, dist);

    if (method_str == "7-4-mer-filter")
        return run2<two_step_k_mer_filter_host_default<7, 4, distance_function>, distance_function>(
            barcodes, read_file_reader, rejection_threshold, dist);

    std::cerr << "Error! Unknown barcode calling method" << method_str << std::endl;
    return 8;
}

int main(int argc, char** argv) {

    // parse arguments
    auto args = parse_arguments(argc, argv);

    if (args.contains("help") || args.contains("h")) {
        print_help();
        return -1;
    }

    barcode_file_str = args.contains("barcodes") ? args["barcodes"] : args["b"];
    read_file_str = args.contains("reads") ? args["reads"] : args["r"];
    threshold_str = get_arg(args, "threshold_distance", "t", "");
    distance_str = get_arg(args, "distance", "d", "sequence-levenshtein");
    verbose = args.count("verbose") || args.count("v");
    gpu = args.count("gpu") || args.count("g");
    method_str = get_arg(args, "method", "m", "4-mer-filter");

    if (barcode_file_str.empty() || read_file_str.empty()) {
        std::cerr << "Error! Missing required arguments --barcodes and --reads" << std::endl;
        return 4;
    }

    // load barcode from disk
    const auto barcodes = barcode_set_bc_reader(barcode_file_str);

    // create a reader object to load chunks of reads from disk
    auto read_file_reader = read_file_fastq_reader(read_file_str);

    if (barcodes.empty()) {
        std::cerr << "Error! Empty barcode file!" << std::endl;
        return 1;
    }

    // read threshold
    int32_t rejection_threshold = DEFAULT_REJECTION_THRESHOLD;
    if (!threshold_str.empty())
        rejection_threshold = std::stoi(threshold_str);

    /****************************************************************************
     * Select the cost and distance measure.
     ***************************************************************************/

    // if no costs are given, use uniform costs
    if (cost_file_str.empty()) {

        if (distance_str == "levenshtein")
            return run1<weighted_levenshtein_v1<unit_costs>>
                (barcodes, read_file_reader, rejection_threshold, weighted_levenshtein_v1<unit_costs>());
        if (distance_str == "sequence-levenshtein")
            return run1<weighted_sequence_levenshtein_v1<unit_costs>>
                (barcodes, read_file_reader, rejection_threshold, weighted_sequence_levenshtein_v1<unit_costs>());

        std::cerr << "Error! Unknown distance!" << std::endl;
        return 3;
    }

    // if explicit costs are given, load the cost file and run with explicit costs
    if (distance_str == "levenshtein")
        return run1<weighted_levenshtein_v1<explicit_costs>>
        (barcodes, read_file_reader, rejection_threshold,
         weighted_levenshtein_v1<explicit_costs>(cost_file_str));
    if (distance_str == "sequence-levenshtein")
        return run1<weighted_sequence_levenshtein_v1<explicit_costs>>
        (barcodes, read_file_reader, rejection_threshold,
         weighted_sequence_levenshtein_v1<explicit_costs>(cost_file_str));

    std::cerr << "Error! Unknown distance!" << std::endl;
    return 3;
}
