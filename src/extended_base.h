//
// Created by agkya on 09.07.25.
//

#ifndef EXTENDED_BASE_H
#define EXTENDED_BASE_H

#include <cassert>
#include <cstdint>

class extended_base {

    /*****************************************************
     * A = 000 = 0
     * C = 001 = 1
     * G = 010 = 2
     * T = 011 = 3
     * _ = 100 = 4
     *****************************************************/

    uint8_t value;

    // dummy constructor
    __host__ __device__ extended_base() : value(UINT8_MAX) {}

public:

    /*extended_base(uint8_t value) :
        value(value) {
        assert(value <= 4);
    }*/

    /**
     * Construct a base from the character 'A', 'C', 'G', 'T', '-'.
     * @param c
     */
    __host__ __device__ extended_base(char c) {

        switch (c) {
        case 'A':
        case 'a':
            value = 0;
            return;
        case 'C':
        case 'c':
            value = 1;
            return;
        case 'G':
        case 'g':
            value = 2;
            return;
        case 'T':
        case 't':
            value = 3;
            return;
        case '_':
        case '-':
            value = 4;
            return;
        default:
            value = UINT8_MAX;
        }
    }

    /**
     * Convert the base to 'A', 'C', 'G', 'T', or '-'.
     */
    __host__ __device__ char to_char() const {
        assert(value <= 4);
        constexpr char table[] = {'A', 'C', 'G', 'T', '-'};
        return table[value];
    }

    /**
     * Convert the base to uint8_t so we can do arithmetics with it.
     */
    __host__ __device__ uint8_t to_uint8() const {
        return value;
    }

};



#endif //EXTENDED_BASE_H
