//
// Created by agkya on 25.02.26.
//

#pragma once

#include "cuda_helper.cuh"
#include "cuda_device.cuh"

namespace quik {

    class read_set_dev {

        size_t* read_start_index_dev = nullptr;
        uint8_t* read_sequence_data_dev = nullptr;
        unsigned* read_lengths_dev = nullptr;
        size_t read_count = 0;

    public:

        /**
         * Initialize the read set as a direct copy of the host reads.
         * @param reads_host
         * @param dev
         * @return
         */
        __host__
        read_set_dev& init(
            const read_set& reads_host,
            const cuda_device& dev = cuda_device(0)) {

            read_count = reads_host.size();

            /************************************************************************************************
             * Allocate gpu memory for the reads
             ************************************************************************************************/

            CUDA_CHECK(cudaSetDevice(dev.id));
            CUDA_CHECK(cudaMalloc(&read_start_index_dev, reads_host.get_start_indices_size()));
            CUDA_CHECK(cudaMalloc(&read_sequence_data_dev, reads_host.get_data_size()));
            CUDA_CHECK(cudaMalloc(&read_lengths_dev, reads_host.get_lengths_size()));

            /************************************************************************************************
            * Copy the reads to gpu memory.
            ************************************************************************************************/

            // read start positions
            CUDA_CHECK(cudaMemcpyAsync(read_start_index_dev, reads_host.get_start_indices(),
                reads_host.get_start_indices_size(), cudaMemcpyHostToDevice, dev.stream));

            // concatenated read lengths
            CUDA_CHECK(cudaMemcpyAsync(read_lengths_dev, reads_host.get_lengths(),
                reads_host.get_lengths_size(), cudaMemcpyHostToDevice, dev.stream));

            // concatenated read sequence data
            CUDA_CHECK(cudaMemcpyAsync(read_sequence_data_dev, reads_host.get_data(),
                reads_host.get_data_size(), cudaMemcpyHostToDevice, dev.stream));

            return *this;
        }


        /**
         * Release all GPU memory.
         */
        __host__
        read_set_dev& finalize() {
            CUDA_CHECK(cudaFree(read_start_index_dev));
            CUDA_CHECK(cudaFree(read_sequence_data_dev));
            CUDA_CHECK(cudaFree(read_lengths_dev));
            return *this;
        }

        __host__ __device__
        size_t size() const {
            return read_count;
        }

        __device__
        unsigned get_read_length(unsigned read_id) const {
            assert(read_id < read_count);
            return read_lengths_dev[read_id];
        }

        __device__
        read operator[](size_t read_id) const {
            assert(read_id < read_count);
            uint8_t* start_ptr = read_sequence_data_dev + read_start_index_dev[read_id];
            unsigned length = read_lengths_dev[read_id];
            /* if (threadIdx.x == 0 && read_id < 10) {
                 printf("read operator[%lu]: read_sequence_data_dev=%p, read_start_index_dev[read_id]=%lu, start_ptr=%p, length=%u, *start_ptr=%u\n"
                     , read_id, read_sequence_data_dev, read_start_index_dev[read_id], start_ptr, length, *start_ptr);
             }*/
            return {start_ptr, length};
        }

        /**
         * Return the base address of the sequence data to the given read.
         * @param read_id
         * @return
         */
        __device__
        const uint8_t* get_read_data_base(size_t read_id) const {
            return &read_sequence_data_dev[read_start_index_dev[read_id]];
        }

        /**
         * Return the first address that does not belong to the read data.
         * @return
         */
        __device__
        const uint8_t* get_read_data_end() const {
            return get_read_data_base(read_count);
        }
    };
}

