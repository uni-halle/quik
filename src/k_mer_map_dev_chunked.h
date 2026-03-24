//
// Created by agkya on 02.03.26.
//

#pragma once

#include "k_mer_map_host.h"
#include "barcode_set.h"
#include "cuda_device.cuh"
#include "cuda_helper.cuh"

namespace quik {
    template<unsigned k, unsigned barcodes_per_chunk>
    class k_mer_map_dev_chunked {

        /*************************************************************************
         * This data structure consists of chunk_count separate k-mer maps,
         * concatenated into two pinned arrays.
         *
         * For each k_mer_map, we have two arrays:
         *
         *    map_entries = unsigned[barcodes.size() * k_mer_start_position_count]
         *    start_indices = size_t[k_mer_count * k_mer_start_position_count+ 1]
         *
         * As we need chunk_count copies, the size of the concatenated arrays are
         * as follows:
         ************************************************************************/

        // host memory
        unsigned *map_entries_host = nullptr;
        size_t *start_indices_host = nullptr;

        // device memory
        unsigned *map_entries_dev = nullptr;
        size_t *start_indices_dev = nullptr;

        static constexpr uint32_t k_mer_count = 1 << 2 * k; // 4^k
        static constexpr size_t k_mer_start_position_count = BARCODE_LENGTH - k + 1;
        static constexpr size_t start_indices_per_chunk = k_mer_count * k_mer_start_position_count;
        static constexpr size_t max_map_entries_per_chunk = barcodes_per_chunk * k_mer_start_position_count;

    public:

        __host__
        k_mer_map_dev_chunked &init(const barcode_set &barcodes, const cuda_device &dev) {
            CUDA_CHECK(cudaSetDevice(dev.id));

            // chunk_count = ceil(barcodes.size() / barcodes_per_chunk)
            const size_t chunk_count = (barcodes.size() + barcodes_per_chunk - 1) / barcodes_per_chunk;

            CUDA_CHECK(cudaMallocHost(&map_entries_host,
                chunk_count * max_map_entries_per_chunk * sizeof(unsigned)));
            CUDA_CHECK(cudaMallocHost(&start_indices_host,
                (chunk_count * start_indices_per_chunk + 1) * sizeof(size_t)));

            /*************************************************************************
             * Partition the barcode set into chunks
             ************************************************************************/

            #pragma omp parallel for
            for (unsigned chunk_id = 0; chunk_id < chunk_count; chunk_id++) {
                size_t barcode_start_id = chunk_id * barcodes_per_chunk;
                size_t barcode_end_id = barcode_start_id + barcodes_per_chunk;

                /*****************************************************************************
                 * Start by computing a k-mer map in host memory that consideres only
                 * the barcodes in the current chunk.
                 ****************************************************************************/

                k_mer_map_host<k> map_chunk(barcodes, barcode_start_id, barcode_end_id);

                /*****************************************************************************
                 * Concatenate the individual lists map_chunk[x][i] into the pinned array
                 *
                 *    map_entries_host[map_entries_per_chunk * chunk_id] to
                 *    map_entries_host[map_entries_per_chunk * (chunk_id+1)]
                 *
                 * The individual list map_chunk[x][i] starts at
                 *
                 *    map_entries_host[map_entries_per_chunk * chunk_id + i * k_mer_count + x]
                 ****************************************************************************/

                size_t running_start_index = max_map_entries_per_chunk * chunk_id;
                start_indices_host[flat(chunk_id, 0, 0)] = running_start_index;

                // for each k-mer start position i within the barcode
                for (unsigned i = 0; i <= BARCODE_LENGTH - k; i++) {
                    // for each possible k-mer x
                    for (unsigned x = 0; x < k_mer_count; x++) {
                        // concat the lists into a single array
                        memcpy(&map_entries_host[running_start_index],
                               map_chunk[x][i].data(),
                               map_chunk[x][i].size() * sizeof(unsigned));
                        running_start_index += map_chunk[x][i].size();
                        start_indices_host[flat(chunk_id, x, i) + 1] = running_start_index;
                    }
                }

                //assert(running_start_index == max_map_entries_per_chunk * (chunk_id+1));
            }

            /***************************************************************************
             * Allocate GPU memory
             **************************************************************************/

            const size_t start_indices_size = (chunk_count * start_indices_per_chunk + 1) * sizeof(size_t);
            const size_t map_entries_size = (chunk_count * max_map_entries_per_chunk) * sizeof(unsigned);
            CUDA_CHECK(cudaMalloc(&map_entries_dev, map_entries_size));
            CUDA_CHECK(cudaMalloc(&start_indices_dev, start_indices_size));

            /***************************************************************************
             * Transfer the k-mer map to the GPU memory
             **************************************************************************/

            CUDA_CHECK(cudaMemcpyAsync(map_entries_dev, map_entries_host,
                map_entries_size, cudaMemcpyHostToDevice, dev.stream));
            CUDA_CHECK(cudaMemcpyAsync(start_indices_dev, start_indices_host,
                start_indices_size, cudaMemcpyHostToDevice, dev.stream));

            return *this;
        }

        k_mer_map_dev_chunked &finalize() {
            CUDA_CHECK(cudaFreeHost(map_entries_host));
            CUDA_CHECK(cudaFreeHost(start_indices_host));
            CUDA_CHECK(cudaFree(map_entries_dev));
            CUDA_CHECK(cudaFree(start_indices_dev));
            return *this;
        }


        /**
         * Return a pointer to the list of barcodes from chunk chunk_id in which
         * the k-mer x starts at position i.
         *
         *
         * @param chunk_id
         * @param k_mer K-mer (binary representation)
         * @param pos Start position from {0,1,..., BARCODE_LENGTH-k-1}
         * @return
         */
        __device__
        const unsigned *get_barcodes(unsigned chunk_id, unsigned k_mer, unsigned pos) const {
            auto start_index = start_indices_dev[flat(chunk_id, k_mer, pos)];
            auto res = &map_entries_dev[start_index];
            return res;
        }

        /**
         * Compress the triple (chunk_id, k_mer, pos) into a single index j such that
         *
         *   the barcode list k_mer_map_host[chunk_id][pos][k_mer] starts at
         *                    k_mer_map_dev_chunks[j].
         *
         * @param chunk_id
         * @param k_mer
         * @param pos
         * @return
         */
        __host__ __device__
        static size_t flat(unsigned chunk_id, unsigned k_mer, unsigned pos) {
            return chunk_id * start_indices_per_chunk + pos * k_mer_count + k_mer;
        }

        /**
         * Return the number of barcodes in which the k-mer x starts at position i.
         *
         * @param chunk_id
         * @param k_mer K-mer (binary representation)
         * @param pos Start position from {0,1,... BARCODE_LENGTH-k-1}
         * @return
         */
        __device__
        size_t get_barcodes_size(unsigned chunk_id, unsigned k_mer, unsigned pos) const {
            auto index = flat(chunk_id, k_mer, pos);
            auto res = start_indices_dev[index + 1] - start_indices_dev[index];
            return res;
        }
    };
}
