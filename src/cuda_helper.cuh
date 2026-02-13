//
// Created by agkya on 03.02.26.
//

#ifndef CUDA_HELPER_CUH
#define CUDA_HELPER_CUH

#include <cuda.h>

#define CUDA_CHECK(x) \
do { \
cudaError_t err = (x); \
if (err != cudaSuccess) { \
fprintf(stderr, "CUDA error %s:%d: %s\n", \
__FILE__, __LINE__, cudaGetErrorString(err)); \
abort(); \
} \
} while (0)

#endif //CUDA_HELPER_CUH
