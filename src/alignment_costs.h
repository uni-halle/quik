//
// Created by agkya on 09.07.25.
//

#ifndef COST_MATRIX_H
#define COST_MATRIX_H

#include <fstream>
#include <sstream>

#include "extended_base.h"

class alignment_costs {

protected:

    int data[5 * 5];

public:
    virtual ~alignment_costs() = default;

    /**
     * Initialize a cost matrix with zero costs for each pair of extended bases.
     */
    alignment_costs() {
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                unsigned ij = i * 5 + j;
                data[ij] = 0;
            }
        }
    }

    alignment_costs(const std::string& cost_file) : alignment_costs() {
        std::ifstream infile(cost_file);
        if (!infile)
            throw std::runtime_error("Could not open cost file " + cost_file);

        char x, y;
        int c_xy;

        while (infile >> x >> y >> c_xy)
            set(x, y, c_xy);
    }

    /**
     * Define the cost of assigning x to y.
     *
     * @param x One of {A,C,G,T,-}.
     * @param y One of {A,C,G,T,-}.
     * @param c_xy cost of assigning x to y.
     * @return
     */
    __host__ __device__ void set(const extended_base x, const extended_base y, const int c_xy) {
        unsigned xy = x.to_uint8() * 5 + y.to_uint8();
        data[xy] = c_xy;
    }

    /**
     * Return the cost of assigning x to y.
     *
     * @param x One of {A,C,G,T,-}.
     * @param y One of {A,C,G,T,-}.
     * @return
     */
    __host__ __device__ int get(const extended_base x, const extended_base y) const {
        unsigned xy = x.to_uint8() * 5 + y.to_uint8();
        return data[xy];
    }

    __host__ virtual std::string to_compact_string() const {
        std::stringstream ss;
        for (char x : {'A', 'C', 'G', 'T', '-'}) {
            for (char y : {'A', 'C', 'G', 'T', '-'}) {
                ss << get(x, y) << ",";
            }
        }
        return ss.str();
    }

};

struct unit_costs : public alignment_costs {
    unit_costs() {
        for (char x : {'A', 'C', 'G', 'T', '-'}) {
            for (char y : {'A', 'C', 'G', 'T', '-'}) {
                if (x != y)
                    set(x, y, 1);
            }
        }
    }

    __host__ virtual std::string to_compact_string() const override {
        std::stringstream ss;
        ss << "<unit-costs>";
        return ss.str();
    }
};

#endif //COST_MATRIX_H
