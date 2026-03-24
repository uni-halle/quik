//
// Created by agkya on 04.03.26.
//

#pragma once

#include "k_mer_filter_host_v2.h"

namespace quik {

    template <unsigned k, typename distance_function>
    using k_mer_filter_host_default = k_mer_filter_host_v2<k, distance_function>;

}
