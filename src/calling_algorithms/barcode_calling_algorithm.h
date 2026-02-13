//
// Created by agkya on 28.01.26.
//

#ifndef BARCODE_CALLING_ALGORITHM_H
#define BARCODE_CALLING_ALGORITHM_H
#include "../barcode_set.h"
#include "../extended_barcode_assignment.h"
#include "../read_set.h"
#include "../distance/distance_measure.h"
#include "../distance/sequence_levenshtein_v4.cuh"

namespace barcode_calling {

    class barcode_calling_algorithm {

    protected:
        const std::string name;
        const distance_measure& dist;

    public:

        barcode_calling_algorithm(const std::string& name,
                                  const distance_measure& dist = sequence_levenshtein_v4())
            : name(name), dist(dist) {};

        virtual ~barcode_calling_algorithm() = default;

        virtual extended_barcode_assignment run(
            const barcode_set& barcodes,
            const read_set& reads) const = 0;

        const std::string& get_name() const {
            return name;
        }
    };

}

#endif //BARCODE_CALLING_ALGORITHM_H
