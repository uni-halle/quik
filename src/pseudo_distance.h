//
// Created by steffen on 16.07.24.
//

#ifndef INC_2OPT_PSEUDO_DISTANCE_H
#define INC_2OPT_PSEUDO_DISTANCE_H

#include "sequence.h"
#include <cassert>
#include <iostream>
#include "constants.h"

template<unsigned k, unsigned window_size>
__host__ __device__ int pseudo_distance(const sequence &b, const sequence &r) {

    int dist = 0;

    /*************************************************************************************************
     * For each position i in the barcode, consider the k-mer x := b[i...i+k-1].
     *
     * For each position j in the interval [i-PSEUDO_DIST_WINDOW_SIZE, ..., i+PSEUDO_DIST_WINDOW_SIZE],
     * consider the k-mer y := r[j...j+k-1].
     *
     * If x matches y, increase the pseudo distance by |i-j|.
     ************************************************************************************************/

    static_assert(k > 0, "invalid value of k");
    static_assert(2 * k < 32, "invalid value of k");
    static constexpr uint32_t mask = (1 << 2 * k) - 1; // bitmask with ones at its 2k least significant positions

    uint32_t x = 0; // bit representation of k-mer x = b[i...i+k-1]
    for (unsigned l = 0; l + 1 < k; l++)
        x = (x << 2) + b[l];

    for (int i = 0; i + k <= sequence::LENGTH; i++) {

        // invariant: x == b[i...i+k-1]
        x = ((x << 2) & mask) + b[i + k - 1];

        int j = 0;
        if (i > window_size)
            j = i - window_size;

        uint32_t y = 0; // bit representation of k-mer y = r[j...j+k-1]
        for (unsigned l = 0; l + 1 < k; l++)
            y = (y << 2) + r[j + l];

        for (; j + k <= sequence::LENGTH && j <= i + window_size; j++) {

            // invariant: y == r[j...j+k-1]
            y = ((y << 2) & mask) + r[j + k - 1];

            // if b[i...i+k-1] == r[j...j+k-1]
            if (x == y) {
                dist += abs(i - j) - sequence::LENGTH;
                //break;
            }
        }
    }

    return dist;
}

int pseudo_distance(const sequence &b, const sequence &r, unsigned k) {
    switch (k) {
        case 1:
            return pseudo_distance<1, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 2:
            return pseudo_distance<2, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 3:
            return pseudo_distance<3, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 4:
            return pseudo_distance<4, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 5:
            return pseudo_distance<5, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 6:
            return pseudo_distance<6, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 7:
            return pseudo_distance<7, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 8:
            return pseudo_distance<8, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 9:
            return pseudo_distance<9, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 10:
            return pseudo_distance<10, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 11:
            return pseudo_distance<11, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 12:
            return pseudo_distance<12, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 13:
            return pseudo_distance<13, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 14:
            return pseudo_distance<14, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        case 15:
            return pseudo_distance<15, PSEUDO_DISTANCE_WINDOW_SIZE>(b, r);
        default:
            std::exit(-2);
    }
}

#endif //INC_2OPT_PSEUDO_DISTANCE_H
