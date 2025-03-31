//
// Created by steffen on 17.05.24.
//

#ifndef INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V5_H
#define INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V5_H

#include "sequence.h"
#include "sequence_levenshtein_distance_v2.cuh"
#include "sequence_levenshtein_distance_v3.cuh"
#include <climits>


/**
 * Calculate min { SL-dist(b,r), upper_bound } which is usually faster than calculating SL-dist(b,r).
 * @param b
 * @param r
 * @param upper_bound
 * @return
 */
__host__ __device__ inline uint8_t
sequence_levenshtein_distance_v5(const sequence &b, const sequence &r, uint8_t upper_bound = UINT8_MAX) {

    /****************************************************************************************************
     * Use the idea by Zorita et al. Starcode: sequence clustering based on all-pairs search.
     * https://academic.oup.com/bioinformatics/article/31/12/1913/213875
     *
     * In the i-th row of the distance matrix. we calculate just the elements in the window
     *    [i-upper_bound, i+upper_bound].
     ***************************************************************************************************/

    if (upper_bound >= sequence::LENGTH)
        return sequence_levenshtein_distance_v3(b, r);

    uint8_t min_dist = upper_bound + 1;   // will be returned

    uint8_t D_i[sequence::LENGTH + 1]; // i-th row of the distance matrix

    /************************************************************************************
     * In row i, we are only interested in the columns from j_start to j_end, where
     *
     *   j_start = max { 0, i - upper_bound - 1} and
     *   j_end   = min { i + upper_bound + 1, sequence::LENGTH }.
     *
     * Thus, we consider a window around the main diagonal D[i][i] which contains at most
     * 2 * upper_bound + 3 cells.
     ************************************************************************************/

    int j_start = 0;

    // initialize the first row (associated to i=0)
    for (int j = j_start; j <= sequence::LENGTH && j <= upper_bound + 1; j++)
        D_i[j] = j;

    // for each other row i=1,2,...sequence::LENGTH
    for (int i = 1; i <= sequence::LENGTH; i++) {

        j_start = max(0, i - upper_bound - 1);
        uint8_t x = D_i[j_start]; // x == D[i-1][j-1]

        // initialize the cells D[i][j_start] and D[i][i+upper_bound+1] (if existing)
        D_i[j_start] = j_start == 0 ? i : upper_bound + 1;
        if (i + upper_bound + 1 <= sequence::LENGTH)
            D_i[i + upper_bound + 1] = upper_bound + 1;

        /*********************************************************************************************
         * For each row i, we calculate the minimum entry min { D[i][j] : j_start <= j <= j_end }.
         * If this minimum exceeds the upper bound, we may return early.
         ********************************************************************************************/

        uint8_t min_i = D_i[j_start];

        // for each cell in between D[i][j_start] and D[i][i+upper_bound+1]
        for (int j = j_start + 1; j <= sequence::LENGTH && j <= i + upper_bound; j++) {

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
             *
             * If D[i][j] is equal to or exceeds the upper bound, we return upper_bound.
             ****************************************************************************/

            D_i[j] = min(y + 1, x + edit_cost);
            D_i[j] = min(D_i[j], D_i[j - 1] + 1);

            min_i = min(min_i, D_i[j]);

            x = y;
        }

        // early termination
        if (min_i >= upper_bound)
            return upper_bound;

        // Sequence-Levenshtein modification
        if (i + upper_bound + 1 >= sequence::LENGTH)  // if D[i][sequence::LENGTH] has been calculated
            min_dist = min(min_dist, D_i[sequence::LENGTH]);
    }

    // Sequence-Levenshtein modification
    for (unsigned j = j_start; j <= sequence::LENGTH; j++)
        min_dist = min(min_dist, D_i[j]);

    return min_dist;
}

#endif //INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V5_H
