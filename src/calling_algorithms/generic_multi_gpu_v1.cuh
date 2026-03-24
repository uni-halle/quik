//
// Created by steffen on 03.07.24.
//

#pragma once

#include "../read.h"
#include "barcode_calling_algorithm.h"
#include "../cuda_device.cuh"
#include "../cuda_helper.cuh"
#include <memory>   // std::unique_ptr

#include "../extended_barcode_assignment_writer.h"

namespace quik {

    template <typename sub_algorithm,
              typename distance_function>
    class generic_multi_gpu_v1 : public barcode_calling_algorithm {

    protected:

        const distance_function distance;

    public:
        generic_multi_gpu_v1(const barcode_set& barcodes,
                             const read_set& reads,
                             const int32_t rejection_threshold = INT32_MAX,
                             const distance_function& dist = distance_function())
            : barcode_calling_algorithm(barcodes, reads, rejection_threshold), distance(dist) {}

        /**
         * Run the algorithm
         **/
        generic_multi_gpu_v1& run() override {

            /*********************************************************************************
             * Determine the number of CUDA capable GPUs.
             *********************************************************************************/

            int device_count = 0;
            CUDA_CHECK(cudaGetDeviceCount(&device_count));

            if (device_count == 0)
                throw std::runtime_error("Error! No GPU found!");

            if (device_count == 1) {

                // create an algorithm on the default device
                sub_algorithm alg(barcodes, reads, rejection_threshold, distance);

                // launch the algorithm and get the result
                alg.run();

                // copy the barcode assignment to the final memory space
                memcpy(closest_barcodes,alg.get_1st_barcodes(),read_count * sizeof(unsigned));
                memcpy(closest_barcodes_2nd,alg.get_2nd_barcodes(),read_count * sizeof(unsigned));
                memcpy(closest_distances,alg.get_1st_distances(),read_count * sizeof(int32_t));
                memcpy(closest_distances_2nd,alg.get_2nd_distances(),read_count * sizeof(int32_t));

                return *this;
            }

            /*********************************************************************************
             * For each device, we have an individual read set and algorithm.
             *********************************************************************************/

            // max_reads_per_device = ceil(reads.size()/device_count)
            const size_t max_reads_per_device = (reads.size() + device_count - 1) / device_count;
            std::vector<read_set> read_chunks;
            read_chunks.reserve(device_count);

            for (int i = 0; i < device_count; i++) {
                size_t read_start_id = i * max_reads_per_device;
                size_t read_end_id = read_start_id + max_reads_per_device;
                if (read_end_id > reads.size())
                    read_end_id = reads.size();
                read_chunks.emplace_back(reads, read_start_id, read_end_id);
            }

            /***********************************************************************
             * Start the algorithm on each GPU
             ************************************************************************/

            std::vector<cuda_device> devices(device_count);
            std::vector<sub_algorithm> sub_algorithms;
            sub_algorithms.reserve(device_count);

            for (int device_id = 0; device_id < device_count; device_id++) {

                // initialize a cuda stream on the device
                devices[device_id].init(device_id);

                // create an algorithm on that device
                sub_algorithms.emplace_back(
                    barcodes, read_chunks[device_id], rejection_threshold,
                    distance, devices[device_id]);

                // launch the calculations (but do not wait for them to finish)
                sub_algorithms[device_id].start();
            }

            /***********************************************************************
             * Wait for the algorithms to finish their work and copy the
             * partial assignments to the right places.
             ************************************************************************/

            for (int device_id = 0; device_id < device_count; device_id++) {

                sub_algorithms[device_id].finalize();

                /*std::cout << "from device " << device_id << std::endl;
                for (unsigned read_id = 0; read_id < 10; read_id++) {
                    std::cout << read_id << "\t"
                     << sub_algorithms[device_id].get_1st_barcodes()[read_id] << "\t"
                    << sub_algorithms[device_id].get_1st_distances()[read_id] << "\t"
                    << sub_algorithms[device_id].get_2nd_barcodes()[read_id] << "\t"
                    << sub_algorithms[device_id].get_2nd_distances()[read_id] << std::endl;
                }*/

                memcpy(closest_barcodes + device_id * max_reads_per_device,
                       sub_algorithms[device_id].get_1st_barcodes(),
                       read_chunks[device_id].size() * sizeof(unsigned));
                memcpy(closest_barcodes_2nd + device_id * max_reads_per_device,
                       sub_algorithms[device_id].get_2nd_barcodes(),
                       read_chunks[device_id].size() * sizeof(unsigned));
                memcpy(closest_distances + device_id * max_reads_per_device,
                       sub_algorithms[device_id].get_1st_distances(),
                       read_chunks[device_id].size() * sizeof(int32_t));
                memcpy(closest_distances_2nd + device_id * max_reads_per_device,
                       sub_algorithms[device_id].get_2nd_distances(),
                       read_chunks[device_id].size() * sizeof(int32_t));

                // free the device
                devices[device_id].finalize();
            }

            return *this;
        }
    };
};
