//
// Created by agkya on 29.01.26.
//

#ifndef BARCODE_ASSIGNMENT_OUTPUT_FORMAT_H
#define BARCODE_ASSIGNMENT_OUTPUT_FORMAT_H

#include <sstream>

#include "barcode_assignment.h"
#include "barcode_set.h"
#include "read_set.h"

namespace barcode_calling {

    class barcode_assignment_writer {

    protected:
        std::stringstream ss;

    public:
        barcode_assignment_writer(
            const barcode_set& barcodes,
            const read_set& reads,
            const barcode_assignment& ass,
            const int rejection_threshold = INT_MAX) {

            //ss << "read\tbarcode\tdistance\n";
            for (unsigned read_id = 0; read_id < reads.size(); ++read_id) {

                // truncate everything right of the first whitespace
                if (ass.get_distance(read_id) <= rejection_threshold) {
                    ss << reads.get_name_of(read_id) << "\t";
                    ss << barcodes.get_name_of(ass.get_assigned_barcode(read_id)) << "\t";
                    ss << ass.get_distance(read_id) << "\n";
                }
            }
        }

        void write(std::ostream& out) const {
            out << ss.str();
        }

    };

}

#endif //BARCODE_ASSIGNMENT_OUTPUT_FORMAT_H
