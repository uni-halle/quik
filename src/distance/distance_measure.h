//
// Created by agkya on 28.01.26.
//

#ifndef DISTANCE_MEASURE_H
#define DISTANCE_MEASURE_H

#include "../barcode.h"
#include "../read.h"

namespace barcode_calling {

    class distance_measure {

    protected:

        const std::string name;

    public:

        __host__ distance_measure(const std::string& name) : name(name) {};

        virtual ~distance_measure() = default;

        __host__ const std::string& get_name() const { return name; }

        __host__ __device__ virtual int32_t operator()(const barcode&, const read&) const = 0;
    };
}

#endif //DISTANCE_MEASURE_H
