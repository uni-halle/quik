//
// Created by steffen on 25.06.24.
//

#pragma once

#include "k_mer_filter_multi_gpu_default.cuh"
#include <string>


namespace quik {

    template <unsigned k_large, unsigned k_small, typename distance_function>
    class two_step_k_mer_filter_multi_gpu_v2 :
    public k_mer_filter_multi_gpu_default<k_large, distance_function> {

        typedef k_mer_filter_multi_gpu_default<k_large, distance_function> super;

    public:

        two_step_k_mer_filter_multi_gpu_v2(
            const barcode_set& barcodes,
            const read_set& reads,
            const int rejection_threshold = INT32_MAX,
            const distance_function& dist = distance_function())
        : super(barcodes, reads, rejection_threshold, dist) {}

        static std::string name() {
            return std::to_string(k_large) + "_" +
                std::to_string(k_small) + "_mer_filter_multi_gpu_v2";
        }

        /**
         * Run the algorithm
         **/
        two_step_k_mer_filter_multi_gpu_v2& run() override {

            /************************************************************************************************
             * Start by calculating a barcode assignment using k_large-mers.
             *
             * For this purpose, just call the super class start method.
             ************************************************************************************************/

            super::run();

            /************************************************************************************************
             * Create a new read set of still unassigned reads.
             ************************************************************************************************/

            std::vector<unsigned> unassigned_read_ids;

            for (unsigned read_id = 0; read_id < super::reads.size(); read_id++) {
                if (super::get_1st_distances()[read_id] > super::rejection_threshold)
                    unassigned_read_ids.push_back(read_id);
            }

            /*******************************************************************************
             * Create a subset of the still unassigned reads.
             ******************************************************************************/

            read_set unassigned_reads_host(super::reads, unassigned_read_ids);

            /************************************************************************************************
             * Use the k-small-mers to assign the unassigned reads.
             ************************************************************************************************/

            k_mer_filter_multi_gpu_default<k_small, distance_function> ass_k_small(
                super::barcodes, unassigned_reads_host, super::rejection_threshold, super::distance);
            ass_k_small.run();

            for (unsigned local_read_id = 0; local_read_id < unassigned_read_ids.size(); local_read_id++) {

                unsigned global_read_id = unassigned_read_ids[local_read_id];

                assert(ass_k_small.is_assigned_to_some_barcode(local_read_id));
                unsigned barcode_id_1 = ass_k_small.get_1st_barcodes()[local_read_id];
                unsigned barcode_id_2 = ass_k_small.get_2nd_barcodes()[local_read_id];
                int32_t distance_1 = ass_k_small.get_1st_distances()[local_read_id];
                int32_t distance_2 = ass_k_small.get_2nd_distances()[local_read_id];

                super::assign_as_1st_barcode(global_read_id, barcode_id_1, distance_1);
                super::assign_as_2nd_barcode(global_read_id, barcode_id_2, distance_2);
            }

            return *this;
        }

    };

}
