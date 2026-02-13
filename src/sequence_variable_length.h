//
// Created by agkya on 27.01.26.
//

#ifndef SEQUENCE_VARIABLE_LENGTH_H
#define SEQUENCE_VARIABLE_LENGTH_H

#include "sequence.h"

class sequence_variable_length : public sequence {

protected:

    unsigned _length = 0;
    uint8_t* _bases = nullptr;

public:

    /**
     * Important: No device constructor
     */
    __host__ sequence_variable_length(const std::string& str)
        : _length(str.length()),
          _bases(new uint8_t[get_raw_data_size()]) {

        // convert four consecutive ASCII characters into a single byte
        for (unsigned i = 0; i < get_raw_data_size(); i++) {

            uint8_t byte = 0;

            // read the four next characters of the string
            unsigned j = 4 * i;
            for (; j < i * 4 + 4 && j < str.size(); j++) {
                byte = (byte << 2) + base(str[j]).to_uint8();
            }
            byte = byte << 2 * (4 * i + 4 - j);
            _bases[i] = byte;
        }
    }

    /**
     *  Create a sequence from pre-allocated sequence data.
     **/
    __device__ sequence_variable_length(uint8_t* bases_data, unsigned length)
        : _length(length), _bases(bases_data) {}

    /**
     * Important: No device destructor!
     */
    __host__~sequence_variable_length() override {
        if (_bases != nullptr)
            delete[] _bases;
    }

    sequence_variable_length(const sequence_variable_length& other) :
        _length(other.length()),
        _bases(new uint8_t[other.get_raw_data_size()]) {
        std::copy(other._bases, other._bases + other.get_raw_data_size(), _bases);
    }

    sequence_variable_length& operator=(const sequence_variable_length& other) {
        if (this != &other) {
            delete[] _bases;
            if (other._bases) {
                _length = other.length();
                _bases = new uint8_t[other.get_raw_data_size()];
                std::copy(other._bases, other._bases + other.get_raw_data_size(), _bases);
            }
        }
        return *this;
    }

    sequence_variable_length(sequence_variable_length&& other) noexcept {
        _bases = other._bases;
        _length = other._length;
        other._bases = nullptr;
        other._length = 0;
    }

    sequence_variable_length& operator=(sequence_variable_length&& other) noexcept {
        if (this != &other) {
            delete[] _bases;

            _bases = other._bases;
            _length = other._length;

            other._bases = nullptr;
            other._length = 0;
        }
        return *this;
    }


    /**
     * Return the binary representation of the i-th nucleotide.
     * @param i
     * @return
     */
    __host__ __device__ base operator[](unsigned i) const override {
        return extract_i_th_nucleotide(_bases, i);
    }

    /**
     * Return the number of bases in this sequence.
     */
    __host__ __device__ unsigned length() const override {
        return _length;
    }

    /**
     * Return the pointer to the sequence data.
     **/
    __host__ const uint8_t* get_raw_data() const {
        return _bases;
    }
};

#endif //SEQUENCE_VARIABLE_LENGTH_H
