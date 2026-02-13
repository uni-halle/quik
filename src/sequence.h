//
// Created by steffen on 17.05.24.
//

#ifndef BARCODE_CALLING_SEQUENCE_H
#define BARCODE_CALLING_SEQUENCE_H

#include <string>
#include <iostream>
#include <cassert>
#include <cstdint>
#include <sstream>
#include "base.h"


class sequence {

public:
    __host__ __device__ virtual ~sequence() = default;
    static constexpr char __nucleotides[] = {'A', 'C', 'G', 'T'};

    static void string_to_bytes(const std::string& str, uint8_t* bases) {

        int array_length = (str.length() + 4 - 1) / 4;

        // convert four consecutive ASCII characters into a single byte
        for (unsigned i = 0; i < array_length; i++) {

            uint8_t byte = 0;

            // read the four next characters of the string
            unsigned j = 4 * i;
            for (; j < i * 4 + 4 && j < str.size(); j++) {
                byte = (byte << 2) + base(str[j]).to_uint8();
            }
            byte = byte << 2 * (4 * i + 4 - j);
            bases[i] = byte;
        }
    }

    /**
     * Return the binary representation of the i-th nucleotide.
     * @param sequence_data Binary representation of sequence data
     * @param index Position between 0 and length()-1
     * @return
     */
    __host__ __device__ static base extract_i_th_nucleotide(const uint8_t* sequence_data, unsigned index) {
        unsigned k = index >> 2; // k = i / 4;
        unsigned r = index - (k << 2); // r = i % 4 such that i = 4k+r
        assert(index == 4*k+r);
        uint8_t c = sequence_data[k];
        uint8_t d = (c >> 2 * (3 - r)) & 3;
        return base::from_uint8(d);
    }

    /**
     * Return the number of bases in this sequence.
     */
    __host__ __device__ virtual unsigned length() const = 0;

    /**
     * Extract the base at a certain position.
     **/
    __host__ __device__ virtual base operator[](unsigned index) const = 0;

    /**
     * Convert to an ASCII string
     */
    __host__ std::string to_string() const {
        std::stringstream ss;
        for (unsigned i = 0; i < length(); i++)
            ss << operator[](i).to_char();
        //str[LENGTH] = '\0';
        return ss.str();
    }

    __host__ __device__ unsigned get_raw_data_size() const {
        // array_length = ceil(length/4) as 4 bases fit into one uint8_t
        unsigned array_length = (length() + 4 - 1) / 4;
        //std::cout << "get_array_length() = " << array_length << std::endl;
        return array_length;
    }
};

#endif //BARCODE_CALLING_SEQUENCE_H
