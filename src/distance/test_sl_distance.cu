//
// Created by steffen on 24.07.24.
//

#include <iostream>
#include "read_sequence_file.h"
#include "read.h"
#include "sequence_levenshtein_distance_v2.cuh"

int main(int argc, char **argv) {

    if (argc != 2) {
        std::cerr << "Usage: ./test_sl_distance <barcode_file>" << std::endl;
        return 0;
    }

    // convert command line arguments
    std::string barcode_file(argv[1]);

    // read barcode file
    auto barcodes = read_sequence_file(barcode_file);

    for (unsigned i = 0; i < barcodes.size(); i++) {
        for (unsigned j = 0; j < barcodes.size(); j++) {

            uint8_t d3 = sequence_levenshtein_distance_v3(barcodes[i], barcodes[j]);
            uint8_t d5 = sequence_levenshtein_distance_v5(barcodes[i], barcodes[j], d3);
            if(d3 != d5) {
                std::cout << "Fehler bei i=" << i << ", j=" << j << ": d3=" << d3 << ", d5=" << d5 << std::endl;
                return -1;
            }
        }
    }

    std::cout << "SUCCESS!" << std::endl;
    return 0;

}
