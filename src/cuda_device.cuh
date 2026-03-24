//
// Created by agkya on 04.03.26.
//

#ifndef BARCODE_CALLING_CUDA_STREAM_CUH
#define BARCODE_CALLING_CUDA_STREAM_CUH
#include "cuda_helper.cuh"

namespace quik {

    struct cuda_device {

        int id = 0;
        cudaStream_t stream = nullptr;

    public:

        cuda_device& init(int device_id) {
            id = device_id;
            CUDA_CHECK(cudaSetDevice(id));
            CUDA_CHECK(cudaStreamCreate(&stream));
            return *this;
        }

        cuda_device& finalize() {
            CUDA_CHECK(cudaStreamDestroy(stream));
            return *this;
        }
    };
}

#endif //BARCODE_CALLING_CUDA_STREAM_CUH
