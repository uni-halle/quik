//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_LEVENSHTEIN_DISTANCE_H
#define BARCODE_CALLING_LEVENSHTEIN_DISTANCE_H

#include "sequence.h"

/**
 * Calculate the Levenshtein distance L-dist(b,r).
 * @param b
 * @param r
 * @return
 */
__host__ __device__ uint8_t
levenshtein_distance_v1(const sequence &b, const sequence &r) {

    uint8_t D_i[sequence::LENGTH + 1]; // i-th row of the distance matrix

    // initialize the first row (associated to i=0)
    for (unsigned j = 0; j <= sequence::LENGTH; j++)
        D_i[j] = j;


    // for each character of the sequences
    for (uint8_t i = 1; i <= sequence::LENGTH; i++) {

        uint8_t x = D_i[0]; // x == D[i-1][j-1]
        D_i[0] = i;
        for (uint8_t j = 1; j <= sequence::LENGTH; j++) {

            uint8_t edit_cost = b[i - 1] == r[j - 1] ? 0 : 1;
            uint8_t y = D_i[j];         // y == D[i-1][j]

            /*****************************************************************************
             * Determine D[i][j] = min {
             *    D[i-1][j-1] + edit_cost,    // replace
             *    D[i-1][j] + 1,              // delete
             *    D[i][j-1] + 1               // insert
             * } = min {
             *    x + edit_cost,              // replace
             *    y + 1,                      // delete
             *    D_i[j-1] + 1                // insert
             * }
             ****************************************************************************/

            D_i[j] = min(y + 1, min(x + edit_cost, D_i[j - 1] + 1));
            x = y;
        }
    }

    return D_i[sequence::LENGTH];
}



#endif //BARCODE_CALLING_LEVENSHTEIN_DISTANCE_H
