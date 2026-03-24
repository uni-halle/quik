//
// Created by agkya on 04.03.26.
//

#pragma once

#include "generic_multi_gpu_v1.cuh"
#include "k_mer_filter_single_gpu_default.cuh"

namespace quik {
    template <unsigned k,
              typename distance_function>
    class k_mer_filter_multi_gpu_v1 : public generic_multi_gpu_v1<
            k_mer_filter_single_gpu_default<k, distance_function>, distance_function> {

        typedef generic_multi_gpu_v1<
            k_mer_filter_single_gpu_default<k, distance_function>, distance_function> base;

    public:

        k_mer_filter_multi_gpu_v1(const barcode_set& barcodes,
                                  const read_set& reads,
                                  const int32_t rejection_threshold = INT32_MAX,
                                  const distance_function& dist = distance_function())
            : base(barcodes, reads, rejection_threshold, dist) {}

        static std::string name() {
            return std::to_string(k) + "_mer_filter_multi_gpu_v1";
        }
    };
}
