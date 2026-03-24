//
// Created by agkya on 29.01.26.
//

#ifndef BARCODE_EXTENDED_ASSIGNMENT_OUTPUT_FORMAT_H
#define BARCODE_EXTENDED_ASSIGNMENT_OUTPUT_FORMAT_H

#include <sstream>

#include "extended_barcode_assignment.h"
#include "barcode_set.h"
#include "read_set.h"

namespace barcode_calling {

    class extended_barcode_assignment_writer {

    protected:
        std::stringstream ss;

    public:
        extended_barcode_assignment_writer(
            const barcode_set& barcodes,
            const read_file& reads,
            const extended_barcode_assignment& ass,
            const int32_t rejection_threshold = INT32_MAX) {

            ss << "read\tbarcode1\tbarcode2\tdistance1\tdistance2\n";
            for (unsigned read_id = 0; read_id < ass.get_read_count(); ++read_id) {

                if (ass.get_1st_distances()[read_id] <= rejection_threshold) {
                    ss << reads.get_fastq_id(read_id) << "\t";
                    //ss << read_id << "\t";
                    //ss << ass.get_1st_barcodes()[read_id] << "\t";
                    //ss << ass.get_2nd_barcodes()[read_id] << "\t";
                    ss << barcodes.get_name_of(ass.get_1st_barcodes()[read_id]) << "\t";
                    ss << barcodes.get_name_of(ass.get_2nd_barcodes()[read_id]) << "\t";
                    ss << ass.get_1st_distances()[read_id] << "\t";
                    ss << ass.get_2nd_distances()[read_id] << "\n";
                }
            }
        }

        void write(std::ostream& out) const {
            out << ss.str();
        }

    };

}

#endif //BARCODE_EXTENDED_ASSIGNMENT_OUTPUT_FORMAT_H
