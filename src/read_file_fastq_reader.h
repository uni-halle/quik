//
// Created by agkya on 29.01.26.
//

#ifndef READ_SET_FASTQ_READER_H
#define READ_SET_FASTQ_READER_H

#include "read_file.h"

namespace barcode_calling {

    class read_file_fastq_reader {

    protected:
        std::ifstream infile;
        size_t processed_read_count = 0;

    public:
        /**
         * Load a read set from a FASTQ file.
         * @param filename
         */
        explicit read_file_fastq_reader(const std::string& filename)
            : infile(filename) {
            if (!infile)
                throw std::runtime_error("Could not open FASTQ file: " + filename);
        }


        /**
         * Load the next chunk of at most read_count reads from the file.
         *
         * @return
         */
        read_file next(size_t read_count = SIZE_MAX) {

            read_file reads;
            std::string line[4];

            for (size_t cnt = 0; cnt < read_count &&
                 std::getline(infile, line[0]) &&
                 std::getline(infile, line[1]) &&
                 std::getline(infile, line[2]) &&
                 std::getline(infile, line[3]); cnt++) {

                // remove @ in line[0]
                line[0].erase(0, 1);

                // truncate everything after the first whitespace
                auto pos = line[0].find_first_of(" \t");
                if (pos != std::string::npos)
                    line[0].erase(pos);

                // remove all characters except 'A,'C','G','T' and their lower case equivalents
                std::erase_if(line[1], [](char c) {
                    return !(c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
                        c == 'a' || c == 'c' || c == 'g' || c == 't');
                });

                reads.add(line[1], line[0]);
                processed_read_count++;
            }

            return reads;
        }

        /**
         * Return the number of reads loaded from the FASTQ file so far.
         * @return
         */
        size_t get_processed_read_count() const {
            return processed_read_count;
        }

    };

}

#endif //READ_SET_FASTQ_READER_H
