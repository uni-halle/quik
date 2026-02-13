//
// Created by steffen on 25.06.24.
//

#ifndef INC_TWO_STEP_KMER_FILTERED_CALLING_GPU_V1_CUH
#define INC_TWO_STEP_KMER_FILTERED_CALLING_GPU_V1_CUH

#include "k_mer_filtered_calling_gpu_v4.cuh"
#include "k_mer_filtered_calling_gpu_v5.cuh"
#include "barcode_calling_algorithm.h"
#include <string>

namespace barcode_calling {

    template <unsigned k_large,
              unsigned k_small>
    class two_step_k_mer_filtered_calling_gpu_v1 : public barcode_calling_algorithm {

    public:

        two_step_k_mer_filtered_calling_gpu_v1(
            const distance_measure& dist,
            unsigned rejection_threshold)
            : barcode_calling_algorithm(std::to_string(k_large) + "_" +
                                        std::to_string(k_small) + "_mer_filtered_calling_gpu_v1",
                                        dist, rejection_threshold) {}

        /**
         * Run the algorithm
         **/
        extended_barcode_assignment run(
            const barcode_set& barcodes,
            const read_set& reads
        ) const override {

            /************************************************************************************************
             * Start by calculating a barcode assignment using k_large-mers
             ************************************************************************************************/

            auto ass_k_large = k_mer_filtered_calling_gpu_v5<k_large>(
                dist, rejection_threshold).run(barcodes, reads);

            /************************************************************************************************
             * Create a new read set of still unassigned reads.
             ************************************************************************************************/

            read_set unassigned_reads;
            std::vector<unsigned> original_read_ids;
            for (unsigned read_id = 0; read_id < reads.size(); read_id++) {
                if (!ass_k_large.is_assigned_to_some_barcode(read_id)) {
                    read r = reads[read_id];
                    original_read_ids.push_back(read_id);
                    unassigned_reads.add(r.to_string(), reads.get_name_of(read_id));
                }
            }

            /************************************************************************************************
             * Use the k-small-mers to assign the unassigned reads.
             ************************************************************************************************/

            auto ass_k_small = k_mer_filtered_calling_gpu_v5<k_small>(
                dist, rejection_threshold).run(barcodes, reads);

            for (unsigned unassigned_read_id = 0; unassigned_read_id < unassigned_reads.size(); unassigned_read_id++) {
                if (ass_k_small.is_assigned_to_some_barcode(unassigned_read_id)) {

                    ass_k_large.assign_as_1st_barcode(unassigned_read_id,
                        ass_k_small.get_1st_barcode(unassigned_read_id),
                        ass_k_small.get_distance_to_1st_barcode(unassigned_read_id)
                        );

                    ass_k_large.assign_as_2nd_barcode(unassigned_read_id,
                        ass_k_small.get_2nd_barcode(unassigned_read_id),
                        ass_k_small.get_distance_to_2nd_barcode(unassigned_read_id)
                        );
                }
            }

            return ass_k_large;
        }

    };

}

#endif //INC_TWO_STEP_KMER_FILTERED_CALLING_GPU_V1_CUH
