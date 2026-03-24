//
// Created by agkya on 04.03.26.
//

#pragma once

#include "k_mer_filter_single_gpu_v09.cuh"

namespace quik {

    template<unsigned k, typename distance_function>
    using k_mer_filter_single_gpu_default = k_mer_filter_single_gpu_v09<k, distance_function>;

}
