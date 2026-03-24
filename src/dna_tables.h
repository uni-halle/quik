//
// Created by agkya on 16.03.26.
//

#ifndef BARCODE_CALLING_DNA_TABLES_H
#define BARCODE_CALLING_DNA_TABLES_H

#include <algorithm>
#include <cstdio>

#include "cuda_helper.cuh"

namespace quik {

    namespace dna_tables {

        /********************************************************************************
         * Translate each value from 0 to 4 to some character:
         *
         *  0 -> A
         *  1 -> C
         *  2 -> T
         *  3 -> G
         *  4 -> -
         *
         * This table is generated at compile time and automatically transferred to
         * each GPU.
         *******************************************************************************/

        static constexpr char uint8_to_char_host[] = {'A', 'C', 'T', 'G', '-'};
        static __device__ __constant__ char uint8_to_char_dev[5] = {'A', 'C', 'T', 'G', '-'};

        /********************************************************************************
         * Translate each ASCII-character to some uint8_t:
         *
         *  A (65), a (97)  -> 0
         *  C (67), c (99)  -> 1
         *  T (84), t (116) -> 2
         *  G (71), g (103) -> 3
         *  - (45), _ (95)  -> 4
         *
         * all other -> UINT8_MAX
         *
         * This table is instantiated at compile time and automatically available from
         * each GPU.
         *******************************************************************************/

        __device__ __constant__ inline uint8_t char_to_uint8_dev[256] = {

#define X UINT8_MAX

         /* 0 1 2 3 4 5 6 7 8 9 */
            X,X,X,X,X,X,X,X,X,X,  // 0-9
            X,X,X,X,X,X,X,X,X,X,  // 10-19
            X,X,X,X,X,X,X,X,X,X,  // 20-29
            X,X,X,X,X,X,X,X,X,X,  // 30-39
            X,X,X,X,X,4,X,X,X,X,  // 40-49
            X,X,X,X,X,X,X,X,X,X,  // 50-59
            X,X,X,X,X,0,X,1,X,X,  // 60-69
            X,3,X,X,X,X,X,X,X,X,  // 70-79
            X,X,X,X,2,X,X,X,X,X,  // 80-89
            X,X,X,X,X,4,X,0,X,1,  // 90-99
            X,X,X,3,X,X,X,X,X,X,  // 100-109
            X,X,X,X,X,X,2,X,X,X,  // 110-119
            X,X,X,X,X,X,X,X,X,X,  // 120-129
            X,X,X,X,X,X,X,X,X,X,  // 130-139
            X,X,X,X,X,X,X,X,X,X,  // 140-149
            X,X,X,X,X,X,X,X,X,X,  // 150-159
            X,X,X,X,X,X,X,X,X,X,  // 160-169
            X,X,X,X,X,X,X,X,X,X,  // 170-179
            X,X,X,X,X,X,X,X,X,X,  // 180-189
            X,X,X,X,X,X,X,X,X,X,  // 190-199
            X,X,X,X,X,X,X,X,X,X,  // 200-209
            X,X,X,X,X,X,X,X,X,X,  // 210-219
            X,X,X,X,X,X,X,X,X,X,  // 220-229
            X,X,X,X,X,X,X,X,X,X,  // 230-239
            X,X,X,X,X,X,X,X,X,X,  // 240-249
            X,X,X,X,X,X           // 250-255

        #undef X
        };

        static constexpr uint8_t char_to_uint8_host[256] = {

#define X UINT8_MAX

            /* 0 1 2 3 4 5 6 7 8 9 */
            X,X,X,X,X,X,X,X,X,X,  // 0-9
            X,X,X,X,X,X,X,X,X,X,  // 10-19
            X,X,X,X,X,X,X,X,X,X,  // 20-29
            X,X,X,X,X,X,X,X,X,X,  // 30-39
            X,X,X,X,X,4,X,X,X,X,  // 40-49
            X,X,X,X,X,X,X,X,X,X,  // 50-59
            X,X,X,X,X,0,X,1,X,X,  // 60-69
            X,3,X,X,X,X,X,X,X,X,  // 70-79
            X,X,X,X,2,X,X,X,X,X,  // 80-89
            X,X,X,X,X,4,X,0,X,1,  // 90-99
            X,X,X,3,X,X,X,X,X,X,  // 100-109
            X,X,X,X,X,X,2,X,X,X,  // 110-119
            X,X,X,X,X,X,X,X,X,X,  // 120-129
            X,X,X,X,X,X,X,X,X,X,  // 130-139
            X,X,X,X,X,X,X,X,X,X,  // 140-149
            X,X,X,X,X,X,X,X,X,X,  // 150-159
            X,X,X,X,X,X,X,X,X,X,  // 160-169
            X,X,X,X,X,X,X,X,X,X,  // 170-179
            X,X,X,X,X,X,X,X,X,X,  // 180-189
            X,X,X,X,X,X,X,X,X,X,  // 190-199
            X,X,X,X,X,X,X,X,X,X,  // 200-209
            X,X,X,X,X,X,X,X,X,X,  // 210-219
            X,X,X,X,X,X,X,X,X,X,  // 220-229
            X,X,X,X,X,X,X,X,X,X,  // 230-239
            X,X,X,X,X,X,X,X,X,X,  // 240-249
            X,X,X,X,X,X           // 250-255

        #undef X
        };

    }

}

#endif //BARCODE_CALLING_DNA_TABLES_H
