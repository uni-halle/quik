//
// Created by steffen on 25.06.24.
//

#ifndef INC_2OPT_BARCODE_H
#define INC_2OPT_BARCODE_H

#include "constants.h"
#include "sequence_fixed_length.h"

class barcode : public sequence_fixed_length<BARCODE_LENGTH> {

public:
    barcode(const std::string& str) : sequence_fixed_length<BARCODE_LENGTH>(str) {}

};


#endif //INC_2OPT_BARCODE_H
