//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_BARCODE_ASSIGNMENT_H
#define BARCODE_CALLING_BARCODE_ASSIGNMENT_H

#include <vector>

#include "cuda_helper.cuh"

namespace quik {

    struct barcode_assignment {

    protected:

        // to each read, we assign one barcode and the distance to the assigned barcode
        unsigned* closest_barcodes = nullptr;
        int32_t* closest_distances = nullptr;

        size_t read_count = 0;

    public:

        /**
         *  Create an empty assignment, in which each read is unassigned.
         */
        barcode_assignment(const size_t read_count) : read_count(read_count) {
            CUDA_CHECK(cudaMallocHost(&closest_barcodes, read_count * sizeof(unsigned)));
            CUDA_CHECK(cudaMallocHost(&closest_distances, read_count * sizeof(int32_t)));
            std::fill_n(closest_barcodes, read_count, UINT_MAX);
            std::fill_n(closest_distances, read_count, INT32_MAX);
        }

        virtual ~barcode_assignment() {
            CUDA_CHECK(cudaFreeHost(closest_barcodes));
            CUDA_CHECK(cudaFreeHost(closest_distances));
        }

        bool is_assigned_to_some_barcode(unsigned read_id) const {
            return closest_barcodes[read_id] != UINT_MAX;
        }

        const unsigned* get_closest_barcodes() const {
            return closest_barcodes;
        }

        const int32_t* get_closest_distances() const {
            return closest_distances;
        }

        void assign_read_to_barcode(unsigned read_id, unsigned barcode_id, int32_t dist) {
            closest_barcodes[read_id] = barcode_id;
            closest_distances[read_id] = dist;
        }

        size_t get_read_count() const {
            return read_count;
        }
    };
}

#endif //BARCODE_CALLING_BARCODE_ASSIGNMENT_H
