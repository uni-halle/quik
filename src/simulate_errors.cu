//
// Created by steffen on 02.07.24.
//

#include <fstream>
#include <iostream>
#include "read.h"
#include "barcode.h"
#include <random>

#include "barcode_assignment.h"
#include "barcode_assignment_writer.h"
#include "barcode_set.h"
#include "barcode_set_bc_reader.h"
#include "read_set_fastq_writer.h"
#include "distance/sequence_levenshtein_v4.cuh"

std::random_device rd;
std::default_random_engine eng(0);
std::uniform_real_distribution<double> distr(0, 1);

using namespace barcode_calling;

// simulate a synthesis error
std::string scramble(const barcode& b,
              double p_ins, double p_del, double p_sub,
              unsigned read_length) {

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
    for (unsigned i = 0; i < b.length(); i++) {

        double p = distr(eng);

        // which type of error occurs?
        if (p < (1 - just_deleted) * p_ins + just_inserted * p_del) { // insertion

            // insert a random base
            ss << nucleotides[(unsigned)(distr(eng) * 4)];
            len++;

            i--; // b[i] is considered again in the next iteration

            just_inserted = 1;
            just_deleted = 0;

        } else if (p < just_deleted * p_ins + (1 - just_inserted) * p_del) { // deletion
            // skip b[i]
            just_inserted = 0;
            just_deleted = 1;
        } else if (p < p_ins + p_del + p_sub) { // substitution

            // sample a random base != b[i]
            unsigned ran;
            do {
                ran = (unsigned)(distr(eng) * 4);
            } while (nucleotides[ran] == b[i].to_char());

            ss << nucleotides[ran];
            len++;

            just_inserted = 0;
            just_deleted = 0;

        } else { // no error
            ss << b[i].to_char();
            len++;

            just_inserted = 0;
            just_deleted = 0;
        }
    }

    // insert random bases until the read length matches its supposed length
    for (; len < read_length; len++)
        ss << nucleotides[(unsigned)(distr(eng) * 4)];

    std::string str = ss.str();

    return str.substr(0, read_length); // truncate additional characters
}

int main(int argc, char** argv) {

    if (argc != 7) {
        std::cerr <<
            "Usage: ./simulate_errors <barcode_file> <prob> <read_count> <read_file> <assignment_file> <read_length>"
            << std::endl;
        std::cerr << "  <barcode_file>" << std::endl;
        std::cerr << "     BC formatted input file containing the barcodes" << std::endl;
        std::cerr << "  <prob>" << std::endl;
        std::cerr << "     Total mutation probability (between 0 and 1)" << std::endl;
        std::cerr << "  <read_count>" << std::endl;
        std::cerr << "     Number of reads to generate from the barcodes" << std::endl;
        std::cerr << "  <read_file>" << std::endl;
        std::cerr << "     Output file to which the read set is written." << std::endl;
        std::cerr << "  <assignment_file>" << std::endl;
        std::cerr << "     Output file to which the ground truth assignment is written" << std::endl;
        std::cerr << "  <read_length>" << std::endl;
        std::cerr << "     Desired length of all reads" << std::endl;
        return 0;
    }

    // convert command line arguments
    std::string barcode_file(argv[1]);
    double p = atof(argv[2]);
    int read_count = atoi(argv[3]);
    std::ofstream read_file(argv[4]);
    std::ofstream label_file(argv[5]);
    unsigned read_length = atoi(argv[6]);

    // read barcode file
    auto barcodes = barcode_set_bc_reader(barcode_file);

    // construct a read set and a ground truth assignment
    read_set reads;
    barcode_assignment ass(read_count);

    // simulate sequencing errors
    for (unsigned read_id = 0; read_id < read_count; ++read_id) {
        unsigned barcode_id = read_id % barcodes.size();
        const barcode& b = barcodes[barcode_id];
        std::string read_sequence = scramble(b, p / 3, p / 3, p / 3, read_length);
        std::string read_name = "read" + std::to_string(read_id + 1);
        read r = reads.add(read_sequence, read_name);
        ass.assign_read_to_barcode(read_id, barcode_id, sequence_levenshtein_v4()(b,r));
    }

    // output read set and ground truth
    read_set_fastq_writer(reads).write(read_file);
    barcode_assignment_writer(barcodes, reads, ass).write(label_file);

    return 0;
}
