//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_RECALL_H
#define BARCODE_CALLING_RECALL_H

#include "barcode_assignment.h"

double recall(const barcode_assignment &assignment, const barcode_assignment &ground_truth) {

    unsigned accepted_reads = 0;

    for (unsigned read_id = 0; read_id < assignment.size(); read_id++) {
        if (assignment[read_id] != UINT_MAX) {
            accepted_reads++;
        }
    }

    return (double) accepted_reads / assignment.size();
}

#endif //BARCODE_CALLING_RECALL_H
