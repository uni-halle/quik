//
// Created by steffen on 17.05.24.
//

#ifndef BARCODE_CALLING_SEQUENCE_H
#define BARCODE_CALLING_SEQUENCE_H

#include <string>
#include <iostream>
#include <cassert>
#include <cstdint>

#include "constants.h"

typedef uint8_t base;

struct sequence {

    static constexpr int LENGTH = SEQUENCE_LENGTH;
    static constexpr int ARRAY_LENGTH =
            (LENGTH + 4 - 1) / 4; // ARRAY_LENGTH = ceil(LENGTH/4) as 4 bases fit into one uint8_t
    base bases[ARRAY_LENGTH];

public:

    /**
     * Dummy constructor.
     */
    __host__ __device__ sequence();

    __host__ sequence(const std::string &seq) {

        assert(seq.size() == LENGTH);

        // convert four consecutive ASCII characters into a single byte
        for (unsigned i = 0; i < ARRAY_LENGTH; i++) {

            uint8_t byte = 0;

            // read the four next characters of the string
            unsigned j = 4 * i;
            for (; j < i * 4 + 4 && j < seq.size(); j++) {
                char c = seq[j];
                uint8_t sign = char_to_base(c);
                byte = (byte << 2) + sign;
            }
            byte = byte << 2 * (4 * i + 4 - j);
            bases[i] = byte;
        }
    }

    /**
     * Return the binary representation of the i-th nucleotide.
     * @param i
     * @return
     */
    __host__ __device__ base operator[](unsigned i) const {

        assert(i < LENGTH);
        unsigned k = i >> 2;        // k = i / 4;
        unsigned r = i - (k << 2);  // r = i % 4 such that i = 4k+r
        assert(i == 4*k+r);
        assert(k < ARRAY_LENGTH);
        uint8_t c = bases[k];
        return (c >> 2 * (3 - r)) & 3;
    }

    __host__ __device__ static base char_to_base(char c) {

        unsigned base[256];
        for (unsigned j = 0; j < 256; j++)
            base[j] = UINT8_MAX; // error value

        base['A'] = 0; // 00
        base['a'] = 0;
        base['C'] = 1; // 01
        base['c'] = 1;
        base['G'] = 2; // 10
        base['g'] = 2;
        base['T'] = 3; // 11
        base['t'] = 3;

        return base[c];
    }

    __host__ __device__ static char base_to_char(uint8_t base) {
        assert(base < 4);
        constexpr char table[] = {'A', 'C', 'G', 'T'};
        return table[base];
    }

    __host__ __device__ void to_string(char* str) const {
        for(unsigned i = 0; i < LENGTH; i++)
            str[i] = base_to_char(operator[](i));
        str[LENGTH] = '\0';
    }
};

inline std::ostream &operator<<(std::ostream &os, const sequence &b) {
    for (unsigned i = 0; i < sequence::LENGTH; i++) {
        os << sequence::base_to_char(b[i]);
    }
    return os;
}

#endif //BARCODE_CALLING_SEQUENCE_H
