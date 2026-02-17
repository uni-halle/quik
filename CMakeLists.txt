cmake_minimum_required(VERSION 3.22)

# --- CUDA Architectures Setup -----------------------------------------

# User override:
#   cmake -DCMAKE_CUDA_ARCHITECTURES=86 ..
#
# Default behavior:
#   - If a GPU is present: build with native architecture
#   - If no GPU is present: fall back to a portable architecture list

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)

    message(STATUS "No CUDA architecture specified by user.")

    # Try native first (works only if a CUDA GPU is visible at configure time)
    set(_TRY_NATIVE TRUE)

    # We detect whether nvidia-smi works (simple heuristic)
    find_program(NVIDIA_SMI nvidia-smi)

    if(NVIDIA_SMI)
        execute_process(
                COMMAND ${NVIDIA_SMI} --query-gpu=compute_cap --format=csv,noheader
                RESULT_VARIABLE _SMI_RESULT
                OUTPUT_VARIABLE _SMI_OUTPUT
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        if(_SMI_RESULT EQUAL 0 AND NOT _SMI_OUTPUT STREQUAL "")
            message(STATUS "CUDA GPU detected (compute capability: ${_SMI_OUTPUT})")
            message(STATUS "Using CMAKE_CUDA_ARCHITECTURES=native")
            set(CMAKE_CUDA_ARCHITECTURES native)
            set(_TRY_NATIVE FALSE)
        endif()
    endif()

    # Fallback if no GPU was detected
    if(_TRY_NATIVE)
        message(WARNING
                "No CUDA GPU detected at configure time.\n"
                "Falling back to a portable architecture set (sm_61, sm_75, sm_86).\n"
                "You can override manually with:\n"
                "  cmake -DCMAKE_CUDA_ARCHITECTURES=XX .."
        )

        # Portable default set
        set(CMAKE_CUDA_ARCHITECTURES "61;75;86")
    endif()

else()
    message(STATUS
            "User provided CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}"
    )
endif()

message(STATUS "Final CUDA archs: ${CMAKE_CUDA_ARCHITECTURES}")

# ----------------------------------------------------------------------

project(quik LANGUAGES CXX CUDA)

set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CXX_STANDARD 20)

# set Release mode as the default
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# configure the barcode length
set(BARCODE_LENGTH 34 CACHE STRING "Sequence length of DNA barcodes")
configure_file(
        ${CMAKE_SOURCE_DIR}/src/constants.h.in
        ${CMAKE_BINARY_DIR}/generated/constants.h
)



find_package(OpenMP REQUIRED)

add_executable(simulate_errors
        src/simulate_errors.cu
)

add_executable(quik
        src/quik.cu
        src/distance/pseudo_distance_all_v1.h
        src/calling_algorithms/k_mer_map_host_v1.h
        src/determine_barcode_length.h
        src/barcode_assignment_writer.h
        src/barcode_set.h
        src/barcode_set_bc_reader.h
        src/cuda_helper.cuh
        src/number_of_assigned_reads.h
)

add_executable(benchmark_barcode_calling
        src/calling_algorithms/benchmark_barcode_calling.cu
        src/determine_barcode_length.h
        src/sequence_fixed_length.h
        src/sequence_variable_length.h
        src/extended_barcode_assignment.h
        src/distance/distance_measure.h
        src/calling_algorithms/barcode_calling_algorithm.h
)

# Aktivieren von OpenMP für den Host-Compiler
if(OpenMP_CXX_FOUND)
    target_link_libraries(quik PUBLIC OpenMP::OpenMP_CXX)
    target_compile_options(quik PUBLIC ${OpenMP_CXX_FLAGS})
    target_link_libraries(benchmark_barcode_calling PUBLIC OpenMP::OpenMP_CXX)
    target_compile_options(benchmark_barcode_calling PUBLIC ${OpenMP_CXX_FLAGS})
endif()

# Aktivieren von OpenMP für den CUDA-Compiler
set_target_properties(quik PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
)

target_compile_options(quik
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>
)

target_compile_options(quik PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--std=c++20>
)

target_compile_options(benchmark_barcode_calling PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--std=c++20>
)

target_compile_options(simulate_errors PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--std=c++20>
)

# use the generated include files with constants

target_include_directories(quik PRIVATE
        ${CMAKE_BINARY_DIR}/generated
)

target_include_directories(benchmark_barcode_calling PRIVATE
        ${CMAKE_BINARY_DIR}/generated
)

target_include_directories(simulate_errors PRIVATE
        ${CMAKE_BINARY_DIR}/generated
)

# relax warning #20013-D
target_compile_options(quik PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
)

target_compile_options(benchmark_barcode_calling PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
)

target_compile_options(simulate_errors PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
)
