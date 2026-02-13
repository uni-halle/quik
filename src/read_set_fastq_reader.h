//
// Created by agkya on 29.01.26.
//

#ifndef READ_SET_FASTQ_READER_H
#define READ_SET_FASTQ_READER_H

#include "read_set.h"

namespace barcode_calling {

    class read_set_fastq_reader : public read_set {

    public:
        /**
         * Load a read set from a FASTQ file.
         * @param filename
         */
        read_set_fastq_reader(const std::string& filename) {

            std::ifstream infile(filename);
            if (!infile)
                throw std::runtime_error("Could not open FASTQ file: " + filename);

            std::string line[4];

            while (std::getline(infile, line[0]) &&
                std::getline(infile, line[1]) &&
                std::getline(infile, line[2]) &&
                std::getline(infile, line[3])) {

                // remove @ in line[0]
                line[0].erase(0,1);

                // truncate everything after the first whitespace
                auto pos = line[0].find_first_of(" \t");
                if (pos != std::string::npos)
                    line[0].erase(pos);

                add(line[1], line[0]);
            }
        }
    };

}

#endif //READ_SET_FASTQ_READER_H
