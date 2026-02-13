//
// Created by agkya on 28.01.26.
//

#ifndef DUMMY_ALGORITHM_WAIT_10MS_H
#define DUMMY_ALGORITHM_WAIT_10MS_H

#include <chrono>
#include <thread>
#include <unistd.h>

#include "barcode_calling_algorithm.h"

class barcode;

namespace barcode_calling {
    class dummy_algorithm_wait_1s : public barcode_calling_algorithm {

    public:
        explicit dummy_algorithm_wait_1s(const distance_measure& dist,
                                         unsigned rejection_threshold)
            : barcode_calling_algorithm("dummy_algorithm_wait_1s", dist, rejection_threshold) {}

        /**
         * Run the algorithm
         **/
        extended_barcode_assignment run(
            const barcode_set& barcodes,
            const read_set& reads
        ) const override {
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            return extended_barcode_assignment(reads.size());
        }
    };
}

#endif //DUMMY_ALGORITHM_WAIT_10MS_H
