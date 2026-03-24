//
// Created by agkya on 09.07.25.
//

#ifndef BASE_H
#define BASE_H

#include "extended_base.h"

namespace quik {

    class base : public extended_base {

        /*****************************************************
         * A = 00 = 0
         * C = 01 = 1
         * T = 10 = 2
         * G = 11 = 3
         *****************************************************/

    public:

        /**
         * Construct a base from the character 'A', 'C', 'G', 'T' and lower case letters.
         *
         * char	ASCII	(c >> 1) & 3
         * A	65	    0
         * C	67	    1
         * G	71    	3
         * T	84	    2
         * a	97	    0
         * c	99	    1
         * g	103	    3
         * t	116	    2
         *
         * @param c
         */
        __host__ __device__
        explicit base(const char c) : extended_base(static_cast<uint8_t>((c >> 1) & 3)) {
            assert(c == 'A' || c == 'C' || c == 'G' || c == 'T'
                || c == 'a' || c == 'c' || c == 'g' || c == 't');
            assert(to_char() == c);
        }

        /**
         * Construct a base from an uint8_t.
         *
         * 0 -> 'A'
         * 1 -> 'C'
         * 2 -> 'T'
         * 3 -> 'G'
         *
         * @param value 0 <= value <= 3
         */
        __host__ __device__
         explicit base(uint8_t value) : extended_base(value) {
            assert(value <= 4);
        }
    };
}

#endif //BASE_H
