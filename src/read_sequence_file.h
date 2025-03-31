//
// Created by steffen on 17.05.24.
//

#ifndef INC_2OPT_READ_SEQUENCE_FILE_H
#define INC_2OPT_READ_SEQUENCE_FILE_H

#include "sequence.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

struct read_sequence_file : public std::vector<sequence> {

    read_sequence_file(const std::string &file) {

        /* Read sequences from file */
        std::ifstream infile(file);
        std::string line;
        while (std::getline(infile, line)) {
            emplace_back(line);
        }

    }

};


#endif //INC_2OPT_READ_SEQUENCE_FILE_H
