//
// Created by agkya on 09.07.25.
//

#ifndef BASE_H
#define BASE_H

#include "extended_base.h"

class base : public extended_base {

    /*****************************************************
     * A = 00 = 0
     * C = 01 = 1
     * G = 10 = 2
     * T = 11 = 3
     *****************************************************/


public:

    /**
     * Construct a base from the character 'A', 'C', 'G', 'T'.
     * @param c
     */
    __host__ __device__ base(const char c) : extended_base(c) {
        assert(c != '-' && c != '_');
    }

    __host__ __device__ static base from_uint8(uint8_t b) {
        static constexpr char table[] = {'A', 'C', 'G', 'T', '-'};
        return base(table[b]);
    }

    __host__ __device__
    bool operator==(const base& b) const {
        return this->to_char() == b.to_char();
    }
};


#endif //BASE_H
