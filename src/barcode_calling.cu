#include <iostream>
#include "read_sequence_file.h"
#include "exact_calling_host.h"
#include "k_mer_filtered_calling_gpu_v1.cuh"

int main(int argc, char **argv) {

    if (argc != 4) {
        std::cout << "Usage: ./kmer_filtered_calling_cpu_v1 <sequence> <reads> <k>" << std::endl;
        return -1;
    }

    std::string barcode_file(argv[1]);
    std::string read_file(argv[2]);
    unsigned k = atoi(argv[3]);

    /* Read sequence and read file */
    std::vector<sequence> barcodes = read_sequence_file(barcode_file);
    std::vector<sequence> reads = read_sequence_file(read_file);

    /* Start sequence calling */
    auto assignment = gpu_filtered_calling_v1(barcodes, reads, k).get_assignment();

    /* Output assignment */
    for (unsigned i = 0; i < reads.size(); i++) {
        if (assignment[i] != (unsigned) -1)
            std::cout << i << " " << assignment[i] << std::endl;
        else
            std::cout << i << " -1" << std::endl;
    }

    return 0;
}
