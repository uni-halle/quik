//
// Created by agkya on 19.03.26.
//

#pragma once

#include <fstream>
#include <iosfwd>
#include <sstream>
#include "extended_base.h"

namespace quik {

    class explicit_costs {

    protected:

        int32_t data[5 * 5] = {
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0
        };

    public:
        /**
         * Initialize a cost matrix with zero costs for each pair of extended bases.
         */
        explicit_costs() = default;

        ~explicit_costs() = default;

        explicit_costs(const std::string& cost_file) {
            std::ifstream infile(cost_file);
            if (!infile)
                throw std::runtime_error("Could not open cost file " + cost_file);

            char x, y;
            int32_t c_xy;

            while (infile >> x >> y >> c_xy)
                set(extended_base(x), extended_base(y), c_xy);
        }

        /**
         * Define the cost of assigning x to y.
         *
         * @param x One of {A,C,G,T,-}.
         * @param y One of {A,C,G,T,-}.
         * @param c_xy cost of assigning x to y.
         * @return
         */
        __host__ __device__
        void set(const extended_base x, const extended_base y, const int32_t c_xy) {
            const unsigned xy = x.to_uint8() * 5 + y.to_uint8();
            assert(xy < 25);
            data[xy] = c_xy;
        }

        /**
         * Return the cost of assigning x to y.
         *
         * @param x One of {A,C,G,T,-}.
         * @param y One of {A,C,G,T,-}.
         * @return
         */
        __host__ __device__
        int32_t get(const extended_base x, const extended_base y) const {
            const unsigned xy = x.to_uint8() * 5 + y.to_uint8();
            return data[xy];
        }

        __host__
        std::string to_string() const {
            std::stringstream ss;
            for (const char x : {'A', 'C', 'G', 'T', '-'}) {
                for (const char y : {'A', 'C', 'G', 'T', '-'}) {
                    ss << get(extended_base(x), extended_base(y)) << ",";
                }
            }
            return ss.str();
        }

        static std::string name() {
            return "explicit_costs";
        }

    };
}
