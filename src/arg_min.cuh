//
// Created by agkya on 24.02.26.
//

#ifndef BARCODE_CALLING_ARG_MIN_CUH
#define BARCODE_CALLING_ARG_MIN_CUH

namespace quik {

    /**
     * Determine from the set of indices an index i such that keys[i] is minimum.
     *
     * This function is supposed to be called from within a CUDA kernel and performs
     * an *in-block-reduction*.
     *
     * @param index_set Subset of indices. (Will be overridden)
     * @param keys Set of keys.
     * @param size Size of index set.
     * @return
     */
    template<unsigned THREADS_PER_BLOCK>
    __device__ unsigned arg_min(unsigned *index_set, const int32_t *keys, unsigned size) {
        unsigned active_elements = size;
        while (active_elements > 1) {
            unsigned half = (active_elements + 1) / 2; // half = ceil(active_elements/2)
            for (unsigned i = threadIdx.x; half + i < active_elements; i += THREADS_PER_BLOCK) {
                unsigned j = half + i;
                if (keys[index_set[i]] > keys[index_set[j]])
                    index_set[i] = index_set[j];
            }
            if (THREADS_PER_BLOCK > 32 )
                __syncthreads();
            active_elements = half;
        }
        return index_set[0];
    }

}

#endif //BARCODE_CALLING_ARG_MIN_CUH