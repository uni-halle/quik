//
// Created by agkya on 04.03.26.
//

#pragma once

#include "generic_multi_gpu_v1.cuh"
#include "k_mer_filter_multi_gpu_v1.cuh"
#include "k_mer_filter_single_gpu_default.cuh"

namespace barcode_calling {

    template <unsigned k,
              typename distance_function>
    using k_mer_filter_multi_gpu_default = k_mer_filter_multi_gpu_v1<k, distance_function>;

}
