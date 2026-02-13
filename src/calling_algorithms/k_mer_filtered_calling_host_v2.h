//
// Created by steffen on 25.06.24.
//

#ifndef INC_2OPT_KMER_FILTERED_CALLING_HOST_V2_H
#define INC_2OPT_KMER_FILTERED_CALLING_HOST_V2_H

#include "../alignment_costs.h"
#include "../extended_barcode_assignment.h"
#include <vector>
#include <chrono>
#include <algorithm>

#include "barcode_calling_algorithm.h"
#include "k_mer_map_host_v1.h"
#include "../distance/distance_measure.h"


namespace barcode_calling {

    template <unsigned k,
              unsigned PSEUDO_DISTANCE_WINDOW_SIZE = 5,
              unsigned MAX_INDEX_SIZE = 100>
    struct k_mer_filtered_calling_host_v2 : public barcode_calling_algorithm {

        explicit k_mer_filtered_calling_host_v2(const distance_measure& dist)
            : barcode_calling_algorithm(std::to_string(k) + "_mer_filtered_calling_host_v2", dist) {}

        /**
         * Run the algorithm
         **/
        extended_barcode_assignment run(
            const barcode_set& barcodes,
            const read_set& reads
        ) const override {

            extended_barcode_assignment ass(reads.size()); // will be returned

            static_assert(2 * k < 32, "invalid value of k");
            static constexpr uint32_t mask = (1 << 2 * k) - 1;
            // bitmask with ones at its 2k least significant positions

            // milliseconds of each phase
            double d[6] = {0, 0, 0, 0, 0, 0};

            /*******************************************************************************************************************
             * Construct a k-mer-map which contains a list of barcode ids for each combination of a
             * k-mer x in {A,C,G,T}^k and each position i in [0...BARCODE_LENGTH-1].
             *
             * Thus, k_mer_map[x][i] is a list of exactly those barcode id's b, for which x == b[i...i+k-1].
             *******************************************************************************************************************/

            auto start = std::chrono::high_resolution_clock::now();

            k_mer_map_host_v1<k> k_mer_map(barcodes);

            d[0] += std::chrono::duration<double, std::milli>(
                std::chrono::high_resolution_clock::now() - start).count();

            /****************************************************************************************************
             * We process reads independent of each other.
             *
             * 1. For each read, we first compute the pseudo distance to *all* barcodes.
             * 2. Then, we determine a subset of barcodes with smallest pseudo distance.
             * 3. For these barcodes, we calculate the SL-distance.
             * 4. Sort the barcodes in the index by their SL-distance.
             ***************************************************************************************************/


#pragma omp parallel for
            for (unsigned read_id = 0; read_id < reads.size(); read_id++) {

                auto my_start = std::chrono::high_resolution_clock::now();

                const read& r = reads[read_id];

                /****************************************************************************************************
                 * Step 1: Compute the pseudo distance to *all* barcodes.
                 ***************************************************************************************************/

                std::vector<int16_t> pseudo_distance(barcodes.size(), 0);
                std::vector<unsigned> barcode_candidates;

                // prepare the first k-1 characters of k-mer m = r[0...k-1]
                uint32_t m = 0; // m = r[0...k-2]
                for (unsigned l = 0; l < k - 1; l++)
                    m = (m << 2) + r[l].to_uint8();

                // for each read position i in which a k-mer may start
                for (int i = 0; i + k <= r.length(); i++) {

                    // invariant: m == r[i...i+k-1]
                    m = ((m << 2) & mask) + r[i + k - 1].to_uint8();

                    // for each position j in a window around position i
                    int j = 0;
                    if (i > PSEUDO_DISTANCE_WINDOW_SIZE)
                        j = i - PSEUDO_DISTANCE_WINDOW_SIZE;
                    for (; j + k <= BARCODE_LENGTH && j <= i + PSEUDO_DISTANCE_WINDOW_SIZE; j++) {

                        // for each barcodes b in which m occurs at position j
                        for (unsigned barcode_id : k_mer_map[m][j]) {
                            if (pseudo_distance[barcode_id] == 0)
                                barcode_candidates.push_back(barcode_id);
                            pseudo_distance[barcode_id] += abs(i - j) - BARCODE_LENGTH;
                        }
                    }
                }

#pragma omp atomic
                d[1] += std::chrono::duration<double, std::milli>(
                    std::chrono::high_resolution_clock::now() - my_start).count();
                my_start = std::chrono::high_resolution_clock::now();

                /****************************************************************************************************
                 * Step 2: Select the barcodes with smallest pseudo distance.
                 ***************************************************************************************************/

                struct index_element {
                    unsigned barcode_id;
                    int16_t distance;
                } index[MAX_INDEX_SIZE];

                unsigned index_size = 0;

                // for each barcode b with non-default pseudo distance q
                for (unsigned barcode_id : barcode_candidates) {
                    int16_t q = pseudo_distance[barcode_id];

                    // insert {b,q} into the sorted index [0...my_index_size-1]
                    unsigned j = index_size;
                    for (; j > 0 && index[j - 1].distance > q; j--) {
                        if (j < MAX_INDEX_SIZE)
                            index[j] = index[j - 1];
                    }

                    if (j < MAX_INDEX_SIZE)
                        index[j] = {barcode_id, q};

                    if (index_size < MAX_INDEX_SIZE)
                        index_size++;
                }

                /*if (read_id == 1) {
                    printf("index of read %i after step 2:\n", read_id);
                    for (unsigned l = 0; l < index_size; l++) {
                        index_element el = index[l];
                        printf("pos=%i, barcode=%i, distance=%i\n", l, el.barcode_id, el.distance);
                    }
                }*/

#pragma omp atomic
                d[2] += std::chrono::duration<double, std::milli>(
                    std::chrono::high_resolution_clock::now() - my_start).count();
                my_start = std::chrono::high_resolution_clock::now();

                /****************************************************************************************************
                 * Step 3: Replace the pseudo distances by the exact distance for each candidate barcode.
                 ***************************************************************************************************/

                for (unsigned i = 0; i < index_size; i++) {
                    const barcode& b = barcodes[index[i].barcode_id];
                    index[i].distance = dist(b, r);
                }

#pragma omp atomic
                d[3] += std::chrono::duration<double, std::milli>(
                    std::chrono::high_resolution_clock::now() - my_start).count();
                my_start = std::chrono::high_resolution_clock::now();

                /****************************************************************************************************
                 * Step 4: Sort the index by the SL-distance
                 ***************************************************************************************************/

                std::sort(index, index + index_size, [](const index_element& x, const index_element& y) {
                    if (x.distance < y.distance)
                        return true;
                    /*if (x.distance == y.distance)
                        return x.barcode_id < y.barcode_id;*/
                    return false;
                });

#pragma omp atomic
                d[4] += std::chrono::duration<double, std::milli>(
                    std::chrono::high_resolution_clock::now() - my_start).count();

                /****************************************************************************************************
                 * Select the barcode with smallest and second-to-smallest SL-distance
                 ***************************************************************************************************/

                ass.assign_as_1st_barcode(read_id, index[0].barcode_id, index[0].distance);
                ass.assign_as_2nd_barcode(read_id, index[1].barcode_id, index[1].distance);
            }

            /****************************************************************************************************
             * Output debug information
             ***************************************************************************************************/

            /* double total = d[0] + d[1] + d[2] + d[3] + d[4];
             std::cout << "Step 0 (preprocessing):       " << 100 * d[0] / total << "%" << std::endl;
             std::cout << "Step 1 (pseudo-distances):    " << 100 * d[1] / total << "%" << std::endl;
             std::cout << "Step 2 (candidate selection): " << 100 * d[2] / total << "%" << std::endl;
             std::cout << "Step 3 (SL-distances):        " << 100 * d[3] / total << "%" << std::endl;
             std::cout << "Step 4 (barcode assignment):  " << 100 * d[4] / total << "%" << std::endl;*/

            return ass;
        }


    };
}

#endif //INC_2OPT_KMER_FILTERED_CALLING_HOST_V2_H
