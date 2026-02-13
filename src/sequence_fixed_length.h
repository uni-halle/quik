//
// Created by agkya on 27.01.26.
//

#ifndef SEQUENCE_FIXED_LENGTH_H
#define SEQUENCE_FIXED_LENGTH_H

#include "sequence.h"

template <unsigned LENGTH>
struct sequence_fixed_length : public sequence {

    // ARRAY_LENGTH = ceil(LENGTH/4) as 4 bases fit into one uint8_t
    static constexpr int ARRAY_LENGTH = (LENGTH + 4 - 1) / 4;
    uint8_t sequence_data[ARRAY_LENGTH];

public:

    sequence_fixed_length() = default;

    __host__ sequence_fixed_length(const std::string& str) {
        assert(str.size() == LENGTH);
        sequence::string_to_bytes(str, sequence_data);
    }

    /**
     * Return the binary representation of the i-th nucleotide.
     * @param i
     * @return
     */
    __host__ __device__ base operator[](unsigned i) const override {
        return extract_i_th_nucleotide(sequence_data, i);
    }

    /**
     * Return the number of bases in this sequence.
     */
    __host__ __device__ unsigned length() const override {
        return LENGTH;
    }
};


#endif //SEQUENCE_FIXED_LENGTH_H
