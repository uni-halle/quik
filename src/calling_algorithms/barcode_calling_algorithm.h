//
// Created by agkya on 28.01.26.
//

#ifndef BARCODE_CALLING_ALGORITHM_H
#define BARCODE_CALLING_ALGORITHM_H

#include "../barcode_set.h"
#include "../extended_barcode_assignment.h"
#include "../read_set.h"

namespace barcode_calling {

    class barcode_calling_algorithm : public extended_barcode_assignment {

    protected:
        const barcode_set& barcodes;
        const read_set& reads;
        const int32_t rejection_threshold = INT32_MAX;

    public:
        /**
         * Initialize the algorithm with the necessary parameters.
         *
         * The actual computation must be started by using the run method.
         *
         * @param barcodes
         * @param reads
         * @param rejection_threshold
         */
        barcode_calling_algorithm(const barcode_set& barcodes,
                                  const read_set& reads,
                                  const int32_t rejection_threshold = INT32_MAX)
            : extended_barcode_assignment(reads.size()),
              barcodes(barcodes), reads(reads), rejection_threshold(rejection_threshold) {};

        /**
         * Run the algorithm.
         * @return
         */
        virtual barcode_calling_algorithm& run() {
            return *this;
        };
    };
}

#endif //BARCODE_CALLING_ALGORITHM_H
