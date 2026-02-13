//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_PRECISION_H
#define BARCODE_CALLING_PRECISION_H

#include "barcode_assignment.h"
#include "read_set.h"

namespace barcode_calling {

    inline double precision(
        const read_set& reads,
        const barcode_assignment& assignment,
        const barcode_assignment& ground_truth,
        const int rejection_threshold = INT_MAX) {

        unsigned accepted_reads = 0;
        unsigned correctly_identified_reads = 0;

        for (unsigned read_id = 0; read_id < reads.size(); read_id++) {
            if (assignment.is_assigned_to_some_barcode(read_id)
                && assignment.get_distance(read_id) <= rejection_threshold) {
                accepted_reads++;
                if (assignment.get_assigned_barcode(read_id) == ground_truth.get_assigned_barcode(read_id))
                    correctly_identified_reads++;
            }
        }

        return static_cast<double>(correctly_identified_reads) / accepted_reads;
    }
}

#endif //BARCODE_CALLING_PRECISION_H
