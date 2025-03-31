//
// Created by steffen on 17.05.24.
//

#ifndef INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V1_H
#define INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V1_H

#include "sequence.h"

__host__ __device__ inline uint8_t sequence_levenshtein_distance_v1(const sequence &b, const sequence &r, const uint8_t upper_bound = UINT8_MAX) {

    uint8_t min_dist = UINT8_MAX;

    uint8_t two_rows[2 * (sequence::LENGTH + 1)]; // two rows of the distance matrix

    uint8_t *current_row = two_rows;
    uint8_t *previous_row = two_rows + sequence::LENGTH + 1;

    // initialize the first row (associated to i=0)
    for (unsigned j = 0; j <= sequence::LENGTH; j++)
        current_row[j] = j;

    // for each character of the sequences
    for (unsigned i = 1; i <= sequence::LENGTH; i++) {

        // swap the roles of current_row and previous_row
        uint8_t *tmp = current_row;
        current_row = previous_row;
        previous_row = tmp;

        current_row[0] = i;
        for (unsigned j = 1; j <= sequence::LENGTH; j++) {
            unsigned edit_cost = b[i - 1] == r[j - 1] ? 0 : 1;

            /*****************************************************************************
             * Determine D[i][j] = min {
             *    D[i-1][j - 1] + edit_cost,  // replace
             *    D[i-1][j] + 1,              // delete
             *    D[i][j - 1] + 1             // insert
             * }
             ****************************************************************************/

            current_row[j] = previous_row[j] + 1; // delete
            if (current_row[j] > previous_row[j - 1] + edit_cost) // replace
                current_row[j] = previous_row[j - 1] + edit_cost;
            if (current_row[j] > current_row[j - 1] + 1) // insert
                current_row[j] = current_row[j - 1] + 1;
        }
        min_dist = min(min_dist, current_row[sequence::LENGTH]);
    }

    for (unsigned j = 0; j <= sequence::LENGTH; j++)
        min_dist = min(min_dist, current_row[j]);

    return min_dist;
}


#endif //INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V1_H
