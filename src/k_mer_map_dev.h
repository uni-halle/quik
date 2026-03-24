//
// Created by agkya on 02.03.26.
//

#ifndef BARCODE_CALLING_K_MER_MAP_DEV_H
#define BARCODE_CALLING_K_MER_MAP_DEV_H

#include "k_mer_map_host.h"
#include "barcode_set.h"
#include "cuda_device.cuh"
#include "cuda_helper.cuh"

namespace barcode_calling {

    template <unsigned k>
    class k_mer_map_dev {

        // host memory
        unsigned* map_entries_host = nullptr;
        size_t* start_indices_host = nullptr;

        // device memory
        unsigned* map_entries_dev = nullptr;
        size_t* start_indices_dev = nullptr;

        static constexpr uint32_t k_mer_count = 1 << 2 * k; // 4^k
        static constexpr unsigned k_mer_start_position_count = BARCODE_LENGTH - k + 1;

    public:
        __host__
        k_mer_map_dev& init(const barcode_set& barcodes,
                            const cuda_device& dev) {

            CUDA_CHECK(cudaSetDevice(dev.id));

            /*****************************************************************************
             * Start by computing the k-mer map in host memory
             ****************************************************************************/

            k_mer_map_host<k> k_mer_map_host(barcodes);

            /*****************************************************************************
             * Create a k-mer map in gpu memory space.
             *
             * For this purpose, we concatenate the individual barcode lists
             * k_mer_map[x][i] into a single array map_entries_dev, similarly to an
             * adjacency array of a graph.
             *
             * The barcode indices from k_map_map[x][i] are then located in the
             * array map_entries_dev at the positions
             *
             *  map_start_index_dev[i * k_mer_count + x] to
             *  map_start_index_dev[i * k_mer_count + x + 1] - 1.
             *
             *  Using cudaMallocHost, we git pinned memory that can be transferred to
             *  the GPU asynchronously.
             ****************************************************************************/

            size_t running_index = 0;

            CUDA_CHECK(cudaMallocHost(&map_entries_host,
                barcodes.size() * k_mer_start_position_count*sizeof(unsigned)));
            CUDA_CHECK(cudaMallocHost(&start_indices_host,
                (k_mer_count * k_mer_start_position_count + 1)*sizeof(size_t)));

            // for each k-mer start position i within the barcode
            for (unsigned i = 0; i <= BARCODE_LENGTH - k; i++) {
                // for each possible k-mer x
                for (unsigned x = 0; x < k_mer_count; x++) {
                    // concat the lists into a single array
                    memcpy(&map_entries_host[running_index], k_mer_map_host[x][i].data(),
                           k_mer_map_host[x][i].size() * sizeof(unsigned));
                    running_index += k_mer_map_host[x][i].size();
                    start_indices_host[i * k_mer_count + x + 1] = running_index;
                }
            }

            /***************************************************************************
             * Transfer the k-mer map to the GPU memory
             **************************************************************************/

            // allocate GPU memory
            constexpr size_t start_indices_size = (k_mer_count * k_mer_start_position_count + 1) * sizeof(size_t);
            const size_t map_entries_size = (barcodes.size() * k_mer_start_position_count) * sizeof(unsigned);
            CUDA_CHECK(cudaMalloc(&map_entries_dev, map_entries_size));
            CUDA_CHECK(cudaMalloc(&start_indices_dev, start_indices_size));

            // copy the data to the GPU memory
            CUDA_CHECK(
                cudaMemcpyAsync(map_entries_dev, map_entries_host,
                    map_entries_size, cudaMemcpyHostToDevice, dev.stream));
            CUDA_CHECK(cudaMemcpyAsync(start_indices_dev, start_indices_host,
                start_indices_size, cudaMemcpyHostToDevice, dev.stream));

            return *this;
        }

        k_mer_map_dev& finalize() {
            CUDA_CHECK(cudaFreeHost(map_entries_host));
            CUDA_CHECK(cudaFreeHost(start_indices_host));
            CUDA_CHECK(cudaFree(map_entries_dev));
            CUDA_CHECK(cudaFree(start_indices_dev));
            return *this;
        }


        /**
         * Return a pointer to the list of barcodes in which the k-mer x starts at position i.
         *
         * This is equivalent to the list k_mer_map[x][i].
         *
         * @param x K-mer (binary representation)
         * @param i Start position from {0,1,..., BARCODE_LENGTH-k-1}
         * @return
         */
        __device__
        const unsigned* get_barcodes(unsigned x, unsigned i) const {
            auto res = &map_entries_dev[start_indices_dev[i * k_mer_count + x]];
            return res;
        }


        /**
         * Return the number of barcodes in which the k-mer x starts at position i.
         * @param x K-mer (binary representation)
         * @param i Start position from {0,1,... BARCODE_LENGTH-k-1}
         * @return
         */
        __device__
        size_t get_barcodes_size(unsigned x, unsigned i) const {
            auto res = start_indices_dev[i * k_mer_count + x + 1] - start_indices_dev[i * k_mer_count + x];
            return res;
        }
    };
}

#endif //BARCODE_CALLING_K_MER_MAP_DEV_H
