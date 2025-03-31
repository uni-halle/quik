//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_BARCODE_ASSIGNMENT_H
#define BARCODE_CALLING_BARCODE_ASSIGNMENT_H

#include <vector>

struct barcode_assignment : public std::vector<unsigned> {

    barcode_assignment(unsigned read_count) :
            std::vector<unsigned>(read_count, UINT_MAX) {};
};

#endif //BARCODE_CALLING_BARCODE_ASSIGNMENT_H
