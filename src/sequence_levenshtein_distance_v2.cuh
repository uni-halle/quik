//
// Created by steffen on 17.05.24.
//

#ifndef INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V2_H
#define INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V2_H

#include "sequence.h"



__host__ __device__ inline uint8_t sequence_levenshtein_distance_v2(const sequence &b, const sequence &r, const uint8_t upper_bound = UINT8_MAX) {

    uint8_t min_dist = UINT8_MAX;

    uint8_t D_i[sequence::LENGTH + 1]; // i-th row of the distance matrix

    // initialize the first row (associated to i=0)
    for (unsigned j = 0; j <= sequence::LENGTH; j++)
        D_i[j] = j;

    // for each character of the sequences
    for (unsigned i = 1; i <= sequence::LENGTH; i++) {

        uint8_t x = D_i[0]; // x == D[i-1][j-1]

        D_i[0] = i;
        for (unsigned j = 1; j <= sequence::LENGTH; j++) {

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

            uint8_t D_ij = y + 1;      // delete
            if (D_ij > x + edit_cost)  // replace
                D_ij = x + edit_cost;
            if (D_ij > D_i[j - 1] + 1) // insert
                D_ij = D_i[j - 1] + 1;
            D_i[j] = D_ij;

            x = y;
        }

        // Sequence-Levenshtein modification
        if (min_dist > D_i[sequence::LENGTH])
            min_dist = D_i[sequence::LENGTH];
    }

    // Sequence-Levenshtein modification
    for (unsigned j = 0; j <= sequence::LENGTH; j++)
        if (min_dist > D_i[j])
            min_dist = D_i[j];

    return min_dist;
}


#endif //INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V2_H
