//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_RECALL_H
#define BARCODE_CALLING_RECALL_H

#include "barcode_assignment.h"

namespace barcode_calling {

    inline double recall(

        const barcode_assignment& assignment,
        unsigned read_count) {

        unsigned accepted_reads = 0;

        for (unsigned read_id = 0; read_id < read_count; read_id++) {
            if (assignment.is_assigned_to_some_barcode(read_id)) {
                accepted_reads++;
            }
        }

        return (double)accepted_reads / read_count;
    }
}

#endif //BARCODE_CALLING_RECALL_H
