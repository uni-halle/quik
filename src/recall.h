//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_RECALL_H
#define BARCODE_CALLING_RECALL_H

#include "barcode_assignment.h"

namespace barcode_calling {

    inline double recall(
        const read_set& reads,
        const barcode_assignment& assignment,
        const int rejection_threshold = INT_MAX) {

        unsigned accepted_reads = 0;

        for (unsigned read_id = 0; read_id < reads.size(); read_id++) {
            if (assignment.is_assigned_to_some_barcode(read_id)
                && assignment.get_distance(read_id) <= rejection_threshold) {
                accepted_reads++;
            }
        }

        return (double)accepted_reads / reads.size();
    }
}

#endif //BARCODE_CALLING_RECALL_H
