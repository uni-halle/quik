//
// Created by agkya on 24.02.26.
//

#ifndef BARCODE_CALLING_BARCODE_ARRAY_H
#define BARCODE_CALLING_BARCODE_ARRAY_H

#include "barcode_set.h"
#include "cuda_device.cuh"

namespace barcode_calling {
    class barcode_set_dev {
        /**
         * Re-alignment of a barcode set in a two-dimensional way:
         *
         *                        barcode_id
         *          |   0  |   1  |   2  |   3  |   4  |  ... | barcodes.size()-1 |
         *          --------------------------------------------------------------|
         *  byte 0  | AACT | AACT |      |      |      |      |                   |
         *  byte 1  | CTTG | CCTG |      |      |      |      |                   |
         *  byte 2  | AACT | AACT |      |      |      |      |                   |
         *  byte 3  | ACCA | TTGC |      |      |      |      |                   |
         *          |                      ....                                   |
         *  byte 8  | CCTG | ACGT |      |      |      |      |                   |
         *          ---------------------------------------------------------------
         *
         * In doing so, threads may read barcodes from global memory efficiently.
         */

        uint8_t *data_host = nullptr;
        uint8_t *data_dev = nullptr;
        size_t pitch = SIZE_MAX;
        size_t barcode_count = 0;

        static constexpr unsigned bytes_per_barcode = (BARCODE_LENGTH + 4 - 1) / 4;

    public:
        /**
         *
         * Re-align a given barcode set into a two-dimensional, memory-efficient way.
         *
         * @param barcodes_host
         * @param dev
         * @return
         */
        __host__
        barcode_set_dev &init(const barcode_set &barcodes_host,
                              const cuda_device &dev = cuda_device(0)) {
            barcode_count = barcodes_host.size();

            // prepare the data layout in host memory
            CUDA_CHECK(cudaMallocHost(&data_host, barcode_count * bytes_per_barcode));
            for (unsigned barcode_id = 0; barcode_id < barcodes_host.size(); ++barcode_id) {
                const barcode &b = barcodes_host[barcode_id];
                for (unsigned i = 0; i < bytes_per_barcode; i++) {
                    data_host[i * barcode_count + barcode_id] = b.data()[i];
                }
            }

            // allocate memory on the device and transfer the host data
            CUDA_CHECK(cudaSetDevice(dev.id));
            CUDA_CHECK(cudaMallocPitch(&data_dev, &pitch, barcode_count, bytes_per_barcode));

            // copy the memory to the gpu
            CUDA_CHECK(cudaMemcpy2DAsync(data_dev ,pitch, data_host,
                barcode_count * sizeof(uint8_t), barcode_count * sizeof(uint8_t),
                bytes_per_barcode, cudaMemcpyHostToDevice, dev.stream));

            return *this;
        }

        __host__
        barcode_set_dev &finalize() {
            CUDA_CHECK(cudaFreeHost(data_host));
            CUDA_CHECK(cudaFree(data_dev));
            return *this;
        }

        __device__ barcode operator[](const unsigned barcode_index) const {
            assert(barcode_index < barcode_count);
            barcode b;
            for (unsigned i = 0; i < bytes_per_barcode; i++) {
                const uint8_t *addr = data_dev + i * pitch + barcode_index * sizeof(uint8_t);
                b.data()[i] = *addr;
            }
            return b;
        }

        /**
         * Return the number of barcodes.
         * @return
         */
        __host__ __device__ size_t size() const {
            return barcode_count;
        }
    };
}

#endif //BARCODE_CALLING_BARCODE_ARRAY_H
