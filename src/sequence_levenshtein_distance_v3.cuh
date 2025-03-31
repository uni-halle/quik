//
// Created by steffen on 17.05.24.
//

#ifndef INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V3_H
#define INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V3_H

#include "sequence.h"

__host__ __device__ inline uint8_t sequence_levenshtein_distance_v3(const sequence &b, const sequence &r, const uint8_t upper_bound = UINT8_MAX) {

    uint8_t min_dist = sequence::LENGTH;
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

            D_i[j] = min(y + 1, x + edit_cost);
            D_i[j] = min(D_i[j], D_i[j - 1] + 1);

            x = y;
        }

        // Sequence-Levenshtein modification
        min_dist = min(min_dist, D_i[sequence::LENGTH]);
    }

    // Sequence-Levenshtein modification
    for (unsigned j = 0; j <= sequence::LENGTH; j++)
        min_dist = min(min_dist, D_i[j]);

    return min_dist;
}


#endif //INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V3_H
