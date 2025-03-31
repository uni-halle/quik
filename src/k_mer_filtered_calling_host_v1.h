//
// Created by steffen on 25.06.24.
//

#ifndef INC_2OPT_KMER_FILTERED_CALLING_HOST_V1_H
#define INC_2OPT_KMER_FILTERED_CALLING_HOST_V1_H

#include "barcode.h"
#include "read.h"
#include "barcode_assignment.h"
#include "index_element.h"
#include "sequence_levenshtein_distance_v2.cuh"
#include <vector>
#include <chrono>
#include <algorithm>

struct k_mer_filtered_calling_host_v1 : public barcode_assignment {

    k_mer_filtered_calling_host_v1(const std::vector<barcode>& barcodes,
                                   const std::vector<read>& reads,
                                   const unsigned k,
                                   const int distance_measure,
                                   const int rejection_threshold = REJECTION_THRESHOLD)
        : barcode_assignment(run(barcodes, reads, k, distance_measure, rejection_threshold)) {}

private:
    barcode_assignment run(const std::vector<barcode>& barcodes,
                           const std::vector<read>& reads,
                           const unsigned k,
                           const int distance_measure,
                           const int rejection_threshold) {
        switch (k) {
        case 4:
            return run<4>(barcodes, reads, distance_measure, rejection_threshold);
        case 5:
            return run<5>(barcodes, reads, distance_measure, rejection_threshold);
        case 6:
            return run<6>(barcodes, reads, distance_measure, rejection_threshold);
        case 7:
            return run<7>(barcodes, reads, distance_measure, rejection_threshold);
        default:
            std::cerr << "invalid parameter k!" << std::endl;
            return barcode_assignment(reads.size()); // return invalid assignment
        }
    }

    template <unsigned k>
    barcode_assignment run(const std::vector<barcode>& barcodes,
                           const std::vector<read>& reads,
                           const int distance_measure,
                           const int rejection_threshold) {
        switch (distance_measure) {
        case SEQUENCE_LEVENSHTEIN_DISTANCE:
            return run<k, sequence_levenshtein_distance_v5>(barcodes, reads, rejection_threshold);
        case LEVENSHTEIN_DISTANCE:
            return run<k, levenshtein_distance_v2>(barcodes, reads, rejection_threshold);
        default:
            std::cerr << "invalid distance measure!" << std::endl;
            return barcode_assignment(reads.size()); // return invalid assignment
        }
    }

    template <
        unsigned k,
        uint8_t (*exact_distance)(const sequence&, const sequence&, uint8_t upper_bound) // distance measure
    >
    barcode_assignment run(const std::vector<barcode>& barcodes,
                                     const std::vector<read>& reads,
                                     const int rejection_threshold) {

        barcode_assignment assigned_barcode(reads.size());

        auto t0 = std::chrono::high_resolution_clock::now();

        /*******************************************************************************************************************
         * Construct a k-mer-map which contains a list of barcode ids for each combination of a
         * k-mer x in {A,C,G,T}^k and each position i in [0...sequence::LENGTH-1].
         *
         * Thus, k_mer_map[x][i] is a list of exactly those barcode id's b, for which x == b[i...i+k-1].
         *******************************************************************************************************************/

        static_assert(2 * k < 32, "Invalid Paramter k");

        static constexpr uint32_t map_size = 1 << 2 * k; // 4^k
        static constexpr uint32_t mask = (1 << 2 * k) - 1; // bitmask with ones at its 2k least significant positions

        std::vector<std::vector<std::vector<unsigned>>> k_mer_map(map_size, std::vector<std::vector<
                                                                      unsigned>>(sequence::LENGTH));

        // for each barcode
        for (unsigned barcode_id = 0; barcode_id < barcodes.size(); barcode_id++) {
            const barcode& b = barcodes[barcode_id];

            uint32_t x = 0; // x = b[i...i+k-1]
            for (unsigned l = 0; l < k - 1; l++)
                x = (x << 2) + b[l];

            for (uint8_t i = 0; i + k <= sequence::LENGTH; i++) {

                // invariant: x == b[i...i+k-1]
                x = ((x << 2) & mask) + b[i + k - 1];

                // k-mer x is located in barcode b at position i
                k_mer_map[x][i].push_back(barcode_id);
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        double d0 = std::chrono::duration<double, std::milli>(t1 - t0).count();
        double d1 = 0;
        double d2 = 0;
        double d3 = 0;
        double d4 = 0;

        /****************************************************************************************************
         * We process reads consecutively.
         *
         * 1. For each read, we first compute the pseudo distance to *all* barcodes.
         * 2. Then, we determine a subset of barcodes with smallest pseudo distance.
         * 3. For these barcodes, we calculate the SL-distance.
         * 4. Sort the barcodes in the index by their SL-distance.
         ***************************************************************************************************/

        for (unsigned read_id = 0; read_id < reads.size(); read_id++) {

            auto t2 = std::chrono::high_resolution_clock::now();

            read r = reads[read_id];

            /****************************************************************************************************
             * Step 1: Compute the pseudo distance to *all* barcodes.
             ***************************************************************************************************/

            unsigned PSEUDO_DISTANCE_DEFAULT = sequence::LENGTH * sequence::LENGTH * sequence::LENGTH;
            std::vector<unsigned> pseudo_distance(barcodes.size(), PSEUDO_DISTANCE_DEFAULT);
            std::vector<unsigned> barcode_candidates;

            uint32_t m = 0; // m = r[i...i+k-1]
            for (unsigned l = 0; l < k - 1; l++)
                m = (m << 2) + r[l];

            // for each k-mer m in the read
            for (int i = 0; i + k <= sequence::LENGTH; i++) {

                // invariant: m == r[i...i+k-1]
                m = ((m << 2) & mask) + r[i + k - 1];

                // for each position j in a window around position i
                int j = 0;
                if (i > PSEUDO_DISTANCE_WINDOW_SIZE)
                    j = i - PSEUDO_DISTANCE_WINDOW_SIZE;
                for (; j + k <= sequence::LENGTH && j <= i + PSEUDO_DISTANCE_WINDOW_SIZE; j++) {

                    // for each barcodes b in which m occurs at position j
                    for (unsigned barcode_id : k_mer_map[m][j]) {
                        if (pseudo_distance[barcode_id] == PSEUDO_DISTANCE_DEFAULT)
                            barcode_candidates.push_back(barcode_id);
                        pseudo_distance[barcode_id] += abs(i - j) - sequence::LENGTH;
                    }
                }
            }

            auto t3 = std::chrono::high_resolution_clock::now();
            d1 += std::chrono::duration<double, std::milli>(t3 - t2).count();

            /****************************************************************************************************
             * Step 2: Select the barcodes with smallest pseudo distance.
             ***************************************************************************************************/

            index_element index[MAX_INDEX_SIZE];
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

            auto t4 = std::chrono::high_resolution_clock::now();
            d2 += std::chrono::duration<double, std::milli>(t4 - t3).count();

            /****************************************************************************************************
             * Step 3: Replace the pseudo distances by SL-distances for each candidate barcode.
             ***************************************************************************************************/

            for (unsigned i = 0; i < index_size; i++) {
                index[i].distance = exact_distance(barcodes[index[i].barcode_id], r, UINT8_MAX);
            }

            auto t5 = std::chrono::high_resolution_clock::now();
            d3 += std::chrono::duration<double, std::milli>(t5 - t4).count();

            /****************************************************************************************************
             * Step 4: Sort the index by the SL-distance
             ***************************************************************************************************/

            std::sort(index, index + index_size, [](const index_element& x, const index_element& y) {
                return x.distance < y.distance;
            });

            auto t6 = std::chrono::high_resolution_clock::now();
            d4 += std::chrono::duration<double, std::milli>(t6 - t5).count();

            /****************************************************************************************************
             * Select the barcode with smallest SL-distance
             ***************************************************************************************************/

            if (index[0].distance <= REJECTION_THRESHOLD)
                assigned_barcode[read_id] = index[0].barcode_id;
        }

        /****************************************************************************************************
         * Output debug information
         ***************************************************************************************************/

        std::chrono::duration<double, std::milli> total = std::chrono::high_resolution_clock::now() - t0;

        /*std::cout << "Step 0: " << d0 << " (" << 100 * d0 / total.count() << "%)" << std::endl;
        std::cout << "Step 1: " << d1 << " (" << 100 * d1 / total.count() << "%)" << std::endl;
        std::cout << "Step 2: " << d2 << " (" << 100 * d2 / total.count() << "%)" << std::endl;
        std::cout << "Step 3: " << d3 << " (" << 100 * d3 / total.count() << "%)" << std::endl;
        std::cout << "Step 4: " << d4 << " (" << 100 * d4 / total.count() << "%)" << std::endl;
        std::cout << "Total : " << total.count() << " (" << 100 * total.count() / total.count() << "%)" << std::endl;*/

        return assigned_barcode;
    }

};

#endif //INC_2OPT_KMER_FILTERED_CALLING_HOST_V1_H
