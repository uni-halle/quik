//
// Created by agkya on 29.01.26.
//

#ifndef BARCODE_SET_FASTQ_READER_H
#define BARCODE_SET_FASTQ_READER_H
#include "barcode_set.h"

namespace barcode_calling {

    class barcode_set_bc_reader : public barcode_set {

    public:
        /**
         * Load a barcode set from a bc file.
         * @param filename
         */
        barcode_set_bc_reader(const std::string& filename) {

            std::ifstream infile(filename);
            if (!infile)
                throw std::runtime_error("Could not open barcode file: " + filename);

            std::string line[2];

            while (std::getline(infile, line[0]) &&
                std::getline(infile, line[1])) {

                // remove @ in line[0]
                line[0].erase(0, 1);

                // truncate everything after the first whitespace
                auto pos = line[0].find_first_of(" \t");
                if (pos != std::string::npos)
                    line[0].erase(pos);

                // test if the second line has BARCODE_LENGTH characters
                if (line[1].length() != BARCODE_LENGTH) {
                    throw std::runtime_error(
                        "Unsupported barcode length " + std::to_string(line[1].length())
                        + ". Use cmake -DBARCODE_LENGTH=" + std::to_string(line[1].length())
                        + " in the build process to recompile the program.");
                }

                add(barcode(line[1]), line[0]);
            }
        }
    };

}

#endif //BARCODE_SET_FASTQ_READER_H
