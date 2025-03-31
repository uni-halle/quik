//
// Created by steffen on 02.07.24.
//

#include <iostream>
#include "read_sequence_file.h"
#include "read.h"
#include "barcode.h"
#include <random>

std::random_device rd;
std::default_random_engine eng(0);
std::uniform_real_distribution<double> distr(0, 1);


// simulate a synthesis error
read scramble(const barcode &b, double p_ins, double p_del, double p_sub) {

    assert(p_ins >= 0 && p_ins <= 1);
    assert(p_del >= 0 && p_del <= 1);
    assert(p_sub >= 0 && p_sub <= 1);
    assert(p_ins + p_del + p_sub <= 1);

    std::stringstream ss;
    unsigned len = 0;

    char nucleotides[] = {'A', 'C', 'G', 'T'};

    int just_inserted = 0;
    int just_deleted = 0;

    // for each position in the barcode
    for (unsigned i = 0; i < barcode::LENGTH; i++) {

        double p = distr(eng);

        // which type of error occurs?
        if (p < (1-just_deleted)*p_ins + just_inserted*p_del) { // insertion

            // insert a random base
            ss << nucleotides[(unsigned) (distr(eng) * 4)];
            len++;

            i--; // b[i] is considered again in the next iteration

            just_inserted = 1;
            just_deleted = 0;

        } else if (p < just_deleted*p_ins + (1-just_inserted)*p_del) { // deletion
            // skip b[i]
            just_inserted = 0;
            just_deleted = 1;
        } else if (p < p_ins + p_del + p_sub) { // substitution

            // sample a random base != b[i]
            unsigned base;
            do {
                base = (unsigned) (distr(eng) * 4);
            } while (base == b[i]);

            ss << nucleotides[base];
            len++;

            just_inserted = 0;
            just_deleted = 0;

        } else { // no error
            ss << nucleotides[b[i]];
            len++;

            just_inserted = 0;
            just_deleted = 0;
        }
    }

    // insert random bases until the read length matches the barcode length
    for (; len < barcode::LENGTH; len++)
        ss << nucleotides[(unsigned) (distr(eng) * 4)];

    std::string str = ss.str();
    assert(str.length() <= len);

    return str.substr(0, sequence::LENGTH); // truncate additional characters
}

int main(int argc, char **argv) {

    if (argc != 6) {
        std::cerr << "Usage: ./simulate_errors <barcode_file> <p> <read_count> <read_file> <label_file>" << std::endl;
        return 0;
    }

    // convert command line arguments
    std::string barcode_file(argv[1]);
    double p = atof(argv[2]);
    int read_count = atoi(argv[3]);
    std::ofstream read_file(argv[4]);
    std::ofstream label_file(argv[5]);

    // read barcode file
    auto barcodes = read_sequence_file(barcode_file);

    // simulate sequencing errors
    while (read_count > 0) {
        for (unsigned barcode_id = 0; barcode_id < barcodes.size() && read_count > 0; barcode_id++) {
            read_file << scramble(barcodes[barcode_id], p / 3, p / 3, p / 3) << std::endl;
            label_file << barcode_id << std::endl;
            read_count--;
        }
    }

    return 0;
}