//
// Created by agkya on 09.07.25.
//

#pragma once

#include <cassert>
#include <cstdint>

#include "dna_tables.h"

namespace barcode_calling {

    class extended_base {

    protected:
        /*****************************************************
         * A = 000 = 0
         * C = 001 = 1
         * T = 010 = 2
         * G = 011 = 3
         * - = 100 = 4
         *****************************************************/

        uint8_t value = UINT8_MAX;

        __host__ __device__
        explicit extended_base(uint8_t value) : value(value) {};

    public:
        /**
         * Construct a base from the character 'A', 'C', 'G', 'T', '-'.
         * @param c
         */
        __host__ __device__ __forceinline__
        explicit extended_base(const char c) : value(
#ifdef __CUDA_ARCH__   // Device Code
            dna_tables::char_to_uint8_dev[c]
#else
            dna_tables::char_to_uint8_host[c]
#endif
        ) {}

        /**
         * Convert the base to 'A', 'C', 'G', 'T', or '-'.
         */
        __host__ __device__ char to_char() const {
            assert(value <= 4);
#ifdef __CUDA_ARCH__   // Device Code
            return dna_tables::uint8_to_char_dev[value];
#else
            return dna_tables::uint8_to_char_host[value];
#endif
        }

        /**
         * Convert the base to uint8_t so we can do arithmetics with it.
         */
        __host__ __device__ uint8_t to_uint8() const {
            return value;
        }

        /*__host__ __device__
        bool operator==(const extended_base& b) const {
            return this->to_char() == b.to_char();
        }*/

    };
}
