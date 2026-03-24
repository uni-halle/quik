//
// Created by steffen on 17.05.24.
//

#pragma once

#include "../barcode.h"
#include "../read.h"
#include "../unit_costs.h"

namespace quik {

    /**
     * Generalized distance which allows for non-uniform costs.
     */
    template <typename alignment_costs = unit_costs>
    class weighted_levenshtein_v1 {

        const alignment_costs c;

    public:
        static std::string name() { return "weighted_levenshtein_v1"; }

        weighted_levenshtein_v1() = default;

        __host__
        weighted_levenshtein_v1(alignment_costs c) : c(c) {}

        const alignment_costs& get_costs() const { return c; }

        __host__ __device__
        int32_t operator()(const barcode& b, const read& r) const {

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
             *    D[i-1][j-1] + c.get(r[i-1], b[j-1])  // substiute b[j-1] by r[i-1]
             *    D[i-1][j] + c.get(r[i-1],'_')        // delete r[i-1]
             *    D[i][j-1] + c.get('_',b[j-1])        // insert b[j-1]
             * }
             *******************************************************************/

            int32_t D_i[BARCODE_LENGTH + 1]; // i-th row of the distance matrix

            // initialize the first row (associated to read position i=0)
            D_i[0] = 0;
            for (unsigned j = 1; j <= BARCODE_LENGTH; j++)
                D_i[j] = D_i[j - 1] + c.get(b[j - 1], extended_base('_'));

            // for each character of the read
            for (unsigned i = 1; i <= r.length(); i++) {

                int32_t x = D_i[0]; // x == D[i-1][j-1]
                D_i[0] = D_i[0] + c.get(extended_base('_'), r[i - 1]);

                for (uint32_t j = 1; j <= BARCODE_LENGTH; j++) {

                    int32_t y = D_i[j]; // y == D[i-1][j]

                    /*****************************************************************************
                     * Determine D[i][j] = min {
                     *    D[i-1][j-1] + c.get(r[i-1], b[j-1])  // substiute r[i-1] by b[j-1]
                     *    D[i-1][j] + c.get(r[i-1],'_')        // delete r[i-1]
                     *    D[i][j-1] + c.get('_',b[j-1])        // insert b[j-1]
                     * } = min {
                     *    x + c.get(r[i-1], b[j-1]),           // substitute r[i-1] by b[j-1]
                     *    y + c.get(r[i-1],'_'),               // delete r[i-1]
                     *    D_i[j-1] + c.get('_',b[j-1])         // insert b[j-1]
                     * }
                     ****************************************************************************/

                    int32_t x1 = x + c.get(r[i - 1], b[j - 1]);
                    int32_t x2 = y + c.get(r[i - 1], extended_base('_'));
                    int32_t x3 = D_i[j - 1] + c.get(extended_base('_'), b[j - 1]);
                    D_i[j] = min(x1, min(x2, x3));

                    x = y;
                }
            }

            return D_i[BARCODE_LENGTH];

        }

    };
}
