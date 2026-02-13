//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_PRECISION_H
#define BARCODE_CALLING_PRECISION_H

#include "barcode_assignment.h"

namespace barcode_calling {

    inline double precision(
        const barcode_assignment& assignment,
        const barcode_assignment& ground_truth,
        unsigned read_count) {

        unsigned accepted_reads = 0;
        unsigned correctly_identified_reads = 0;

        for (unsigned read_id = 0; read_id < read_count; read_id++) {
            if (assignment.is_assigned_to_some_barcode(read_id)) {
                accepted_reads++;
                if (assignment.get_assigned_barcode(read_id) == ground_truth.get_assigned_barcode(read_id))
                    correctly_identified_reads++;
            }
        }

        return static_cast<double>(correctly_identified_reads) / accepted_reads;
    }
}

#endif //BARCODE_CALLING_PRECISION_H
