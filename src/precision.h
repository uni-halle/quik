//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_PRECISION_H
#define BARCODE_CALLING_PRECISION_H

#include "barcode_assignment.h"

double precision(const barcode_assignment &assignment, const barcode_assignment &ground_truth) {

    unsigned accepted_reads = 0;
    unsigned correctly_identified_reads = 0;

    for (unsigned read_id = 0; read_id < assignment.size(); read_id++) {
        if (assignment[read_id] != UINT_MAX) {
            accepted_reads++;
            if (assignment[read_id] == ground_truth[read_id])
                correctly_identified_reads++;
        }
    }

    return (double) correctly_identified_reads / accepted_reads;
}

#endif //BARCODE_CALLING_PRECISION_H
