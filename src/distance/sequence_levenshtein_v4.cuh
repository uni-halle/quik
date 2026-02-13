//
// Created by steffen on 17.05.24.
//

#ifndef INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V4_H
#define INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V4_H

#include "../barcode.h"
#include "../read.h"
#include "sequence_levenshtein_v3.cuh"

namespace barcode_calling {


    /**
    * Calculate min { SL-dist(b,r), upper_bound } which is usually faster than calculating SL-dist(b,r).
    * @return
    */
    class sequence_levenshtein_v4 : public distance_measure {

    protected:
        uint8_t upper_bound = UINT8_MAX;

    public:
        __host__ __device__ static uint8_t evaluate(const barcode& b, const read& r, uint8_t upper_bound = UINT8_MAX) {

            /************************************************************************************
             * We construct a matrix D with r.length() + 1 rows and BARCODE_LENGTH + 1 columns.
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
             ************************************************************************************/

            assert((int) BARCODE_LENGTH <= (int) UINT8_MAX);
            assert((int) r.length() <= (int) UINT8_MAX);

            /****************************************************************************************************
             * Use the idea by Zorita et al. Starcode: sequence clustering based on all-pairs search.
             * https://academic.oup.com/bioinformatics/article/31/12/1913/213875
             *
             * In the i-th row of the distance matrix. we calculate just the elements in the window
             *    [i-upper_bound, i+upper_bound].
             ***************************************************************************************************/

            // fall back to the standard DP-algorithm
            if (upper_bound >= b.length() || upper_bound >= r.length())
                return sequence_levenshtein_v3::evaluate(b, r);

            uint8_t min_dist = upper_bound; // will be returned
            uint8_t D_i[BARCODE_LENGTH + 1]; // i-th row of the distance matrix

            /************************************************************************************
             * In row i, we are only interested in the columns from j_start to j_end, where
             *
             *   j_start = max { 0, i - upper_bound - 1} and
             *   j_end   = min { i + upper_bound + 1, BARCODE_LENGTH }.
             *
             * Thus, we consider a window around the main diagonal D[i][i] which contains at most
             * 2 * upper_bound + 3 cells.
             ************************************************************************************/

            int j_start = 0;
            int j_end = min(upper_bound + 1, BARCODE_LENGTH);

            // initialize the first row (associated to i=0)
            for (int j = j_start; j <= j_end; j++)
                D_i[j] = j;

            // for each other row i=1,2,...r.length()
            for (int i = 1; i <= r.length(); i++) {

                j_start = max(0, i - upper_bound - 1);
                uint8_t x = D_i[j_start]; // x == D[i-1][j-1]

                // initialize the cells D[i][j_start] and D[i][i+upper_bound+1] (if existing)
                D_i[j_start] = j_start == 0 ? i : upper_bound + 1;
                if (i + upper_bound + 1 <= BARCODE_LENGTH)
                    D_i[i + upper_bound + 1] = upper_bound + 1;

                // for each cell in between D[i][j_start] and D[i][i+upper_bound+1]
                for (int j = j_start + 1; j <= BARCODE_LENGTH && j <= i + upper_bound; j++) {

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
                if (i + upper_bound + 1 >= BARCODE_LENGTH) // if D[i][BARCODE_LENGTH] has been calculated
                    min_dist = min(min_dist, D_i[BARCODE_LENGTH]);
            }

            // Sequence-Levenshtein modification
            for (unsigned j = j_start; j <= BARCODE_LENGTH; j++)
                min_dist = min(min_dist, D_i[j]);

            return min_dist;
        }

        sequence_levenshtein_v4() : distance_measure("sequence_levenshtein_v4") {}

        __host__ __device__
        void set_upper_bound(uint8_t upper_bound) { this->upper_bound = upper_bound; }

        __host__ __device__ int32_t
        operator()(const barcode& b, const read& r) const override {
            return evaluate(b, r, upper_bound);
        }

    };
}

#endif //INC_2OPT_SEQUENCE_LEVENSHTEIN_DISTANCE_V4_H
