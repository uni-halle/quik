//
// Created by steffen on 25.06.24.
//

#ifndef INC_2OPT_READ_H
#define INC_2OPT_READ_H

#include "sequence.h"

namespace barcode_calling {

    class read : public sequence {

    public:

        const uint8_t* sequence_data = nullptr;
        unsigned sequence_length = UINT_MAX;

    public:

        /**
         * Create a read from pre-allocated sequence data.
         *
         * @param sequence_data byte-coded sequence data
         * @param length Number of ASCII characters.
         */
        __host__ __device__
        read(const uint8_t* sequence_data, unsigned length)
            : sequence_data(sequence_data), sequence_length(length) {}

        __host__ __device__
        read& operator=(const read& r) {
            sequence_data = r.sequence_data;
            sequence_length = r.sequence_length;
            return *this;
        }

        __host__ __device__
        unsigned length() const override {
            return sequence_length;
        };

        __host__ __device__
        base operator[](unsigned index) const override {
            return extract_i_th_nucleotide(sequence_data, index);
        }

        __host__ __device__
        const uint8_t* data() const {
            return sequence_data;
        }

    };

}

#endif //INC_2OPT_READ_H
