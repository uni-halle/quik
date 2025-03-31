//
// Created by steffen on 25.06.24.
//

#ifndef INC_2OPT_KMER_FILTERED_CALLING_HOST_V3_CUH
#define INC_2OPT_KMER_FILTERED_CALLING_HOST_V3_CUH

#include "barcode.h"
#include "read.h"
#include "barcode_assignment.h"
#include "k_mer_filtered_calling_host_v2.h"
#include <vector>


struct two_step_k_mer_filtered_calling_host_v1 : public barcode_assignment {

    two_step_k_mer_filtered_calling_host_v1(
            const std::vector <barcode> &barcodes,
            const std::vector <read> &reads,
            const unsigned k_small,
            const unsigned k_large,
            const int distance_measure,
            const int rejection_threshold_small,
            const int rejection_threshold_large)
            : barcode_assignment(run(barcodes, reads, k_small, k_large, distance_measure,rejection_threshold_small,rejection_threshold_large)) {}

private:

    barcode_assignment run(const std::vector <barcode> &barcodes,
                           const std::vector <read> &reads,
                           const unsigned k_small,
                           const unsigned k_large,
                           const int distance_measure,
                           const int rejection_threshold_small,
                           const int rejection_threshold_large
                           ) {

        /************************************************************************************************
         * Start by calculating a barcode assignment using 7-mer-filter.
         ************************************************************************************************/

        barcode_assignment ass_k_large = k_mer_filtered_calling_host_v2(
            barcodes, reads, k_large, distance_measure, rejection_threshold_large);

        /************************************************************************************************
         * Create a list of still unassigned reads.
         ************************************************************************************************/

        std::vector<read> unassigned_reads;
        std::vector<unsigned> original_read_ids;
        for (unsigned read_id = 0; read_id < reads.size(); read_id++) {
            if(ass_k_large[read_id] == UINT_MAX) {
                original_read_ids.push_back(read_id);
                unassigned_reads.push_back(reads[read_id]);
            }
        }

        /************************************************************************************************
         * Use the k-mer filtered method to assign the unassigned reads.
         ************************************************************************************************/

        barcode_assignment ass_k_small = k_mer_filtered_calling_host_v2(
            barcodes, unassigned_reads, k_small, distance_measure, rejection_threshold_small);

        for(unsigned unassigned_read_id = 0; unassigned_read_id < unassigned_reads.size(); unassigned_read_id++) {
            if(ass_k_small[unassigned_read_id] != UINT_MAX) {
                ass_k_large[original_read_ids[unassigned_read_id]] = ass_k_small[unassigned_read_id];
            }
        }

        return ass_k_large;
    }

};

#endif //INC_2OPT_KMER_FILTERED_CALLING_HOST_V3_CUH
