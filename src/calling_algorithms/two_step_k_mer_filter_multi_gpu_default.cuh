//
// Created by agkya on 04.03.26.
//

#pragma once
#include "two_step_k_mer_filter_multi_gpu_v1.cuh"

namespace barcode_calling {

    template <unsigned k_large, unsigned k_small, typename distance_function>
    using two_step_k_mer_filter_multi_gpu_default = two_step_k_mer_filter_multi_gpu_v1<
        k_large, k_small, distance_function>;


}
