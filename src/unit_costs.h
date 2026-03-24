//
// Created by agkya on 09.07.25.
//

#pragma once

#include <fstream>

#include "base.h"

namespace quik {

    struct unit_costs {

        unit_costs() = default;

        ~unit_costs() = default;

        /**
         * Return the cost of assigning x to y.
         *
         * @param x One of {A,C,G,T,-}.
         * @param y One of {A,C,G,T,-}.
         * @return
         */
        __host__ __device__ __forceinline__
        int32_t get(const extended_base x, const extended_base y) const {
            assert(x.to_uint8()*y.to_uint8() < 24); // forbid (-,-)
            return (x.to_uint8() == y.to_uint8()) ? 0 : 1;
        }

        __host__
        std::string to_string() const {
            return "unit_costs";
        }

        static std::string name() {
            return "unit_costs";
        }

    };
}
