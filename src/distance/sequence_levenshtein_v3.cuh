//
// Created by steffen on 17.05.24.
//

#ifndef INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V3_H
#define INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V3_H

#include "distance_measure.h"

namespace barcode_calling {

    /**
     * Calculate the SL-distance between some read and barcode.
     */
    class sequence_levenshtein_v3 : public distance_measure {

    public:

        __host__ __device__ static uint8_t evaluate(const barcode& b, const read& r) {

            /********************************************************************
             * We construct a matrix D with r.length() + 1 rows and
             * BARCODE_LENGTH + 1 columns.
             *
             * The i-th row is associated to the i-1-th base r[i-1] of the read.
             * The j-th column is associated to the j-1-th base b[j-1] of the barcode.
             *
             * The matrix is defined by the following recursion:
             *
             * D[i][j] = min {
             *    D[i-1][j-1] + [0 or 1]  // substitute b[j-1] by r[i-1]
             *    D[i-1][j] + 1           // delete r[i-1]
             *    D[i][j-1] + 1           // insert b[j-1]
             * }
             *
             * We determine the minimum entry of the matrix' last row and column.
             *
             * By definition, the entries D[i][i] are zero for all i.
             * Consequently, no entry is larger than max{r.length(), BARCODE_LENGTH}.
             * We assume than r and b are small enough such that
             *     max{r.length(), BARCODE_LENGTH} <= UINT8_MAX.
             *******************************************************************/

            assert((int) BARCODE_LENGTH <= (int) UINT8_MAX);
            assert((int) r.length() <= (int) UINT8_MAX);

            uint8_t min_dist = UINT8_MAX; // minimum entry of distance matrix D
            uint8_t D_i[BARCODE_LENGTH + 1]; // i-th row of the distance matrix

            // initialize the first row (associated to i=0)
            for (unsigned j = 0; j <= BARCODE_LENGTH; j++)
                D_i[j] = j;

            // for each character of the sequences
            for (uint8_t i = 1; i <= r.length(); i++) {

                uint8_t x = D_i[0]; // x == D[i-1][j-1]
                D_i[0] = i;

                for (uint8_t j = 1; j <= BARCODE_LENGTH; j++) {

                    uint8_t edit_cost = b[i - 1] == r[j - 1] ? 0 : 1;
                    uint8_t y = D_i[j]; // y == D[i-1][j]

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
                min_dist = min(min_dist, D_i[BARCODE_LENGTH]);
            }

            // Sequence-Levenshtein modification
            for (unsigned j = 0; j <= BARCODE_LENGTH; j++)
                min_dist = min(min_dist, D_i[j]);

            return min_dist;
        }

    public:

        sequence_levenshtein_v3() : distance_measure("sequence_levenshtein_v3") {}

        __host__ __device__ inline int32_t operator()(const barcode& b, const read& r) const override {
            return evaluate(b, r);
        }

    };

}

#endif //INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V3_H
