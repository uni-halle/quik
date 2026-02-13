//
// Created by agkya on 13.02.26.
//

#ifndef NUMBER_OF_ASSIGNED_READS_H
#define NUMBER_OF_ASSIGNED_READS_H

#include "read_set.h"
#include "barcode_assignment.h"

namespace barcode_calling {

    inline unsigned number_of_assigned_reads(
        const read_set& reads,
        const barcode_assignment& ass,
        const int rejection_threshold) {
        unsigned num = 0;
        for (unsigned read_id = 0; read_id < reads.size(); ++read_id) {
            if (ass.is_assigned_to_some_barcode(read_id) && ass.get_distance(read_id) <= rejection_threshold)
                num++;
        }
        return num;
    }

}

#endif //NUMBER_OF_ASSIGNED_READS_H
