//
// Created by agkya on 28.01.26.
//

#ifndef EXTENDED_BARCODE_ASSIGNMENT_H
#define EXTENDED_BARCODE_ASSIGNMENT_H

#include "barcode_assignment.h"
#include <cuda.h>

#include "cuda_helper.cuh"

namespace barcode_calling {

    class extended_barcode_assignment : public barcode_assignment {

    protected:

        /**
         * To each read r, we assign two different barcodes:
         *
         *   1) A barcode b1 to which r has minimum distance.
         *   2) A barcode b2 to which r has second best distance.
         */
        unsigned* closest_barcodes_2nd = nullptr;
        int32_t* closest_distances_2nd = nullptr;

    public:

        extended_barcode_assignment(const unsigned read_count)
            : barcode_assignment(read_count) {
            CUDA_CHECK(cudaMallocHost(&closest_barcodes_2nd, read_count * sizeof(unsigned)));
            CUDA_CHECK(cudaMallocHost(&closest_distances_2nd, read_count * sizeof(int32_t)));
            std::fill_n(closest_barcodes_2nd, read_count, UINT_MAX);
            std::fill_n(closest_distances_2nd, read_count, INT32_MAX);
        }

        ~extended_barcode_assignment() override {
            CUDA_CHECK(cudaFreeHost(closest_barcodes_2nd));
            CUDA_CHECK(cudaFreeHost(closest_distances_2nd));
        }

        /**
         * Define the best-fitting barcode to some read.
         * @param read_id
         * @param barcode_id
         * @param distance
         */
        void assign_as_1st_barcode(unsigned read_id, unsigned barcode_id, int distance) {
            assign_read_to_barcode(read_id, barcode_id, distance);
        }

        /**
         * Assign the second-best-fitting barcode to some read.
         * @param read_id
         * @param barcode_id
         * @param distance
         */
        void assign_as_2nd_barcode(unsigned read_id, unsigned barcode_id, int distance) {
            closest_barcodes_2nd[read_id] = barcode_id;
            closest_distances_2nd[read_id] = distance;;
        }


        /**
        * Return an barcode index minimum distance to the read.
        * @return
        */
        [[nodiscard]]
        const unsigned* get_1st_barcodes() const {
            return get_closest_barcodes();
        }

        /**
         * Return an barcode index second-to-minimum distance to the read.
         * @return
         */
        [[nodiscard]]
        const unsigned* get_2nd_barcodes() const {
            return closest_barcodes_2nd;
        }

        [[nodiscard]]
        const int32_t* get_1st_distances() const {
            return closest_distances;
        }

        [[nodiscard]]
        const int32_t* get_2nd_distances() const {
            return closest_distances_2nd;
        }

    };
}

#endif //EXTENDED_BARCODE_ASSIGNMENT_H
