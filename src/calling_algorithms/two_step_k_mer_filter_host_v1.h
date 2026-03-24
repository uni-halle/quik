//
// Created by steffen on 25.06.24.
//

#ifndef INC_TWO_STEP_KMER_FILTERED_CALLING_HOST_V1_CUH
#define INC_TWO_STEP_KMER_FILTERED_CALLING_HOST_V1_CUH

#include "k_mer_filter_host_default.cuh"
#include "barcode_calling_algorithm.h"
#include <string>

namespace quik {

    template <unsigned k_large,
              unsigned k_small,
              typename distance_function>
    class two_step_k_mer_filter_host_v1 : public barcode_calling_algorithm {

        const distance_function distance;

    public:
        two_step_k_mer_filter_host_v1(
            const barcode_set& barcodes,
            const read_set& reads,
            const int32_t rejection_threshold = INT32_MAX,
            const distance_function& dist = distance_function())
            : barcode_calling_algorithm(barcodes, reads, rejection_threshold),
              distance(dist) {}

        static std::string name() {
            return std::to_string(k_large) + "_" +
                std::to_string(k_small) + "_mer_filter_host_v1";
        }

        /**
         * Run the algorithm
         **/
        two_step_k_mer_filter_host_v1& run() override {

            /************************************************************************************************
             * Start by calculating a barcode assignment using k_large-mers
             ************************************************************************************************/

            k_mer_filter_host_default<k_large, distance_function> ass_k_large(
                barcodes, reads, rejection_threshold, distance);
            ass_k_large.run();

            /************************************************************************************************
             * Create a new read set of still unassigned reads.
             ************************************************************************************************/

            read_file unassigned_read_file;
            std::vector<unsigned> original_read_ids;

            for (unsigned read_id = 0; read_id < reads.size(); read_id++) {

                read r = reads[read_id];
                assert(ass_k_large.is_assigned_to_some_barcode(read_id));

                unsigned barcode_id_1 = ass_k_large.get_1st_barcodes()[read_id];
                unsigned barcode_id_2 = ass_k_large.get_2nd_barcodes()[read_id];
                int32_t dist_1 = ass_k_large.get_1st_distances()[read_id];
                int32_t dist_2 = ass_k_large.get_2nd_distances()[read_id];

                if (dist_1 > rejection_threshold) {
                    original_read_ids.push_back(read_id);
                    unassigned_read_file.add(r.to_string());
                } else {
                    assign_as_1st_barcode(read_id, barcode_id_1, dist_1);
                    assign_as_2nd_barcode(read_id, barcode_id_2, dist_2);
                }
            }

            // convert the file of unassigned reads into a read set
            read_set unassigned_reads(unassigned_read_file);

            /************************************************************************************************
             * Use the k-small-mers to assign the unassigned reads.
             ************************************************************************************************/

            k_mer_filter_host_default<k_small, distance_function> ass_k_small(
                barcodes, unassigned_reads, rejection_threshold, distance);
            ass_k_small.run();

            for (unsigned unassigned_read_id = 0;
                 unassigned_read_id < unassigned_reads.size();
                 unassigned_read_id++) {

                assert(ass_k_small.is_assigned_to_some_barcode(unassigned_read_id));
                unsigned barcode_id_1 = ass_k_small.get_1st_barcodes()[unassigned_read_id];
                unsigned barcode_id_2 = ass_k_small.get_2nd_barcodes()[unassigned_read_id];
                int32_t distance_1 = ass_k_small.get_1st_distances()[unassigned_read_id];
                int32_t distance_2 = ass_k_small.get_2nd_distances()[unassigned_read_id];

                unsigned read_id = original_read_ids[unassigned_read_id];
                assign_as_1st_barcode(read_id, barcode_id_1, distance_1);
                assign_as_2nd_barcode(read_id, barcode_id_2, distance_2);
            }

            return *this;
        }

    };

}

#endif //INC_TWO_STEP_KMER_FILTERED_CALLING_HOST_V1_CUH
