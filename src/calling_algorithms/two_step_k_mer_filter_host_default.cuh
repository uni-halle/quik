//
// Created by agkya on 04.03.26.
//

#pragma once

#include "two_step_k_mer_filter_host_v1.h"

namespace quik {

    template <unsigned k_large,
              unsigned k_small,
              typename distance_function>
    using two_step_k_mer_filter_host_default = two_step_k_mer_filter_host_v1<k_large, k_small, distance_function>;

}
