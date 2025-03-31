//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_LEVENSHTEIN_DISTANCE_V3_H
#define BARCODE_CALLING_LEVENSHTEIN_DISTANCE_V3_H

#include "sequence.h"
#include "levenshtein_distance_v1.cuh"

/**
 * Calculate min { L-dist(b,r), upper_bound } which is usually faster than calculating L-dist(b,r).
 * @param b
 * @param r
 * @return
 */
__host__ __device__ uint8_t
levenshtein_distance_v2(const sequence &b, const sequence &r, uint8_t upper_bound = UINT8_MAX) {

    /****************************************************************************************************
     * Use the idea by Zorita et al. Starcode: sequence clustering based on all-pairs search.
     * https://academic.oup.com/bioinformatics/article/31/12/1913/213875
     *
     * In the i-th row of the distance matrix. we calculate just the elements in the window
     *    [i-upper_bound, i+upper_bound].
     ***************************************************************************************************/

    if (upper_bound >= sequence::LENGTH)
        return levenshtein_distance_v1(b, r);

    uint8_t D_i[sequence::LENGTH + 1]; // i-th row of the distance matrix

    // initialize the first row (associated to i=0)
    for (uint8_t j = 0; j <= sequence::LENGTH && j <= upper_bound + 1; j++)
        D_i[j] = j;

    /**********************************************************************+
     * Output
     **********************************************************************/

    /*printf("Levenshtein distance v2 (upper_bound=%u):\n", upper_bound);
    printf("       ");
    for (unsigned j = 0; j < sequence::LENGTH; j++)
        printf(" %2c", sequence::base_to_char(r[j]));
    printf("\n    ");
    for (uint8_t j = 0; j <= sequence::LENGTH && j <= upper_bound + 1; j++)
        printf(" %2u", D_i[j]);
    for (uint8_t j = upper_bound + 2; j <= sequence::LENGTH; j++)   // empty spaces to fill the row
        printf("   ");*/

    // for each character of the sequences
    for (uint8_t i = 1; i <= sequence::LENGTH; i++) {

        uint8_t j_start = 0;
        uint8_t x = D_i[j_start]; // x == D[i-1][j-1]

        // initialization of row i
        D_i[j_start] = i;
        if (i > upper_bound + 1) {
            j_start = i - (upper_bound + 1);
            D_i[j_start] = upper_bound + 1;
        }
        if (i + upper_bound + 1 <= sequence::LENGTH)
            D_i[i + upper_bound + 1] = upper_bound + 1;

        for (unsigned j = j_start + 1; j <= sequence::LENGTH && j <= i + upper_bound; j++) {

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

            D_i[j] = min(y + 1, min(x + edit_cost, D_i[j - 1] + 1));

            if (D_i[j] >= upper_bound)
                return upper_bound;

            x = y;
        }

        /***********************************************************************************
         * Output
         ***********************************************************************************/

        /*printf("\n %2c ", sequence::base_to_char(b[i - 1]));
        for (uint8_t j = 0; j < j_start; j++)   // empty spaces to fill the row
            printf("   ");
        for (uint8_t j = j_start; j <= i + upper_bound + 1 && j <= sequence::LENGTH; j++)
            printf(" %2u", D_i[j]);
        for (uint8_t j = i + upper_bound + 2; j <= sequence::LENGTH; j++)   // empty spaces to fill the row
            printf("   ");*/

    }
    //printf("\n");

    return D_i[sequence::LENGTH];
}


#endif //BARCODE_CALLING_LEVENSHTEIN_DISTANCE_V3_H
