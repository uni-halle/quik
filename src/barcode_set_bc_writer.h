//
// Created by agkya on 29.01.26.
//

#ifndef BARCODE_SET_FASTQ_WRITER_H
#define BARCODE_SET_FASTQ_WRITER_H

#include "barcode_set.h"
#include <sstream>
#include <ostream>

namespace barcode_calling {

    class barcode_set_fastq_writer {

      protected:

        std::stringstream ss;

        /**
         * Convert a barcode set to the bc format.
.        *
         * @param filename
         */
        barcode_set_fastq_writer(const barcode_set& barcodes) {
          for(unsigned barcode_id = 0; barcode_id < barcodes.size(); ++barcode_id) {
            ss << barcodes.get_identifier_of(barcode_id)<< '\n';
            ss << barcodes[barcode_id] << '\n';
          }
        }

        void write(std::ostream& out) const {
            out << ss.str();
        }

    };

}

#endif //BARCODE_SET_FASTQ_WRITER_H
