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

struct read_label_file : public barcode_assignment {

    read_label_file(const std::string &file, unsigned read_cnt)
      : barcode_assignment(read_cnt) {

        /* Read labels from file */
        std::ifstream infile(file);
        std::string line;
        for(unsigned read_id = 0; read_id < read_cnt; ++read_id) {

          if(!std::getline(infile, line))
            std::cerr << "Error! corrupted label file" << std::endl;

          int barcode_id = std::stoi(line);
          operator[](read_id) = barcode_id;
        }
    }

};


#endif //INC_2OPT_READ_LABEL_FILE_H
