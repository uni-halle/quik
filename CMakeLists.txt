cmake_minimum_required(VERSION 3.22)
project(quik LANGUAGES CXX CUDA)

set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_ARCHITECTURES 80)

find_package(OpenMP REQUIRED)

add_executable(simulate_errors
        src/simulate_errors.cu
)

add_executable(benchmark_barcode_calling
        src/benchmark_barcode_calling.cu
        src/two_step_k_mer_filtered_calling_gpu_v1.cuh
)

add_executable(benchmark_sl_distance
        src/benchmark_sl_distance.cu
)

# Aktivieren von OpenMP für den Host-Compiler
if(OpenMP_CXX_FOUND)
    target_link_libraries(benchmark_barcode_calling PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(benchmark_sl_distance PUBLIC OpenMP::OpenMP_CXX)
    target_compile_options(benchmark_barcode_calling PUBLIC ${OpenMP_CXX_FLAGS})
    target_compile_options(benchmark_sl_distance PUBLIC ${OpenMP_CXX_FLAGS})
endif()

# Aktivieren von OpenMP für den CUDA-Compiler
set_target_properties(benchmark_barcode_calling PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
)

# Aktivieren von OpenMP für den CUDA-Compiler
set_target_properties(benchmark_sl_distance PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
)

target_compile_options(benchmark_barcode_calling
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>
)

target_compile_options(benchmark_sl_distance
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>
)
