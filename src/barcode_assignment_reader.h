//
// Created by steffen on 17.05.24.
//

#ifndef INC_2OPT_READ_LABEL_FILE_H
#define INC_2OPT_READ_LABEL_FILE_H

#include "barcode_assignment.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace barcode_calling {

    class barcode_assignment_reader : public barcode_assignment {

    public:

        barcode_assignment_reader(const std::string& filename,
                                  const barcode_set& barcodes,
                                  const read_set& reads)
            : barcode_assignment(reads.size()) {

            std::ifstream infile(filename);
            std::string line;

            // read header
            if (!std::getline(infile, line))
                throw std::runtime_error("Error! corrupted assignment file");

            /* Read labels from file */
            for (unsigned i = 0; i < reads.size(); ++i) {

                if (!std::getline(infile, line))
                    throw std::runtime_error("Error! corrupted assignment file");

                auto parts = split_semicolon_or_tab(line);
                if (parts.size() != 3)
                    throw std::runtime_error("Error! corrupted assignment file");

                std::string read_name = parts[0];
                std::string barcode_name = parts[1];
                int dist = std::stoi(parts[2]);

                // determine the index of read and barcode
                unsigned read_id = reads.get_index_of(parts[0]);
                unsigned barcode_id = barcodes.get_index_of(parts[1]);
                assign_read_to_barcode(read_id, barcode_id, dist);
            }
        }

    private:

        std::vector<std::string> split_semicolon_or_tab(const std::string& s) {
            std::vector<std::string> result;

            std::size_t start = 0;
            while (true) {
                std::size_t pos = s.find_first_of(";\t", start);

                if (pos == std::string::npos) {
                    result.emplace_back(s.substr(start));
                    break;
                }

                result.emplace_back(s.substr(start, pos - start));
                start = pos + 1;
            }

            return result;
        }


    };

}
#endif //INC_2OPT_READ_LABEL_FILE_H
