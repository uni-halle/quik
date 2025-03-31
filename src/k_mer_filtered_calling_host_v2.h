//
// Created by steffen on 25.06.24.
//

#ifndef INC_2OPT_KMER_FILTERED_CALLING_HOST_V2_H
#define INC_2OPT_KMER_FILTERED_CALLING_HOST_V2_H

#include "barcode.h"
#include "read.h"
#include "barcode_assignment.h"
#include "sequence_levenshtein_distance_v2.cuh"
#include <vector>
#include <chrono>
#include <algorithm>

struct k_mer_filtered_calling_host_v2 : public barcode_assignment {

public:
    k_mer_filtered_calling_host_v2(const std::vector<barcode>& barcodes,
                                   const std::vector<read>& reads,
                                   const unsigned k,
                                   const int distance_measure,
                                   const int rejection_threshold = REJECTION_THRESHOLD
    )
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
            std::cerr << "Unsupported value k=" << k << std::endl;
            return barcode_assignment(reads.size());
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

        /*******************************************************************************************************************
         * Construct a k-mer-map which contains a list of barcode ids for each combination of a
         * k-mer x in {A,C,G,T}^k and each position i in [0...sequence::LENGTH-1].
         *
         * Thus, k_mer_map[x][i] is a list of exactly those barcode id's b, for which x == b[i...i+k-1].
         *******************************************************************************************************************/

        static_assert(2 * k < 32, "invalid value of k");
        static constexpr uint32_t map_size = 1 << 2 * k; // 4^k
        static constexpr uint32_t mask = (1 << 2 * k) - 1; // bitmask with ones at its 2k least significant positions

        double d0 = 0;
        double d1 = 0;
        double d2 = 0;
        double d3 = 0;
        double d4 = 0;
        auto t0 = std::chrono::high_resolution_clock::now();

        // create empty k-mer-map
        std::vector<std::vector<std::vector<unsigned>>> k_mer_map(map_size,
                                                                  std::vector<std::vector<unsigned>>(sequence::LENGTH));

        // for each k-mer position i
#pragma omp parallel for
        for (unsigned i = 0; i <= sequence::LENGTH - k; i++) {

            // for each barcode
            for (unsigned barcode_id = 0; barcode_id < barcodes.size(); barcode_id++) {
                const barcode& b = barcodes[barcode_id];

                // consider the k-mer x = b[i...i+k-1]
                uint32_t x = 0; // x = b[i...i+k-1]
                for (unsigned l = 0; l < k; l++)
                    x = (x << 2) + b[i + l];
                assert(x < map_size);

                // k-mer x is located in barcode b at position i
                k_mer_map[x][i].push_back(barcode_id);
            }
        }

        /****************************************************************************************************
         * We process reads consecutively.
         *
         * 1. For each read, we first compute the pseudo distance to *all* barcodes.
         * 2. Then, we determine a subset of barcodes with smallest pseudo distance.
         * 3. For these barcodes, we calculate the SL-distance.
         * 4. Sort the barcodes in the index by their SL-distance.
         ***************************************************************************************************/

        auto t1 = std::chrono::high_resolution_clock::now();
        d0 += std::chrono::duration<double, std::milli>(t1 - t0).count();

#pragma omp parallel for
        for (unsigned read_id = 0; read_id < reads.size(); read_id++) {

            auto t0 = std::chrono::high_resolution_clock::now();

            const read& r = reads[read_id];

            /****************************************************************************************************
             * Step 1: Compute the pseudo distance to *all* barcodes.
             ***************************************************************************************************/

            std::vector<int16_t> pseudo_distance(barcodes.size(), 0);
            std::vector<unsigned> barcode_candidates;

            // prepare the first k-1 characters of k-mer m = r[0...k-1]
            uint32_t m = 0; // m = r[0...k-2]
            for (unsigned l = 0; l < k - 1; l++)
                m = (m << 2) + r[l];

            // for each position i in which a k-mer may start
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
                        if (pseudo_distance[barcode_id] == 0)
                            barcode_candidates.push_back(barcode_id);
                        pseudo_distance[barcode_id] += abs(i - j) - sequence::LENGTH;
                    }
                }
            }

#pragma omp atomic
            d1 += std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();

            /****************************************************************************************************
             * Step 2: Select the barcodes with smallest pseudo distance.
             ***************************************************************************************************/

            t0 = std::chrono::high_resolution_clock::now();

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
            d2 += std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();

            /****************************************************************************************************
             * Step 3: Replace the pseudo distances by SL-distances for each candidate barcode.
             ***************************************************************************************************/

            t0 = std::chrono::high_resolution_clock::now();

            for (unsigned i = 0; i < index_size; i++) {
                const barcode& b = barcodes[index[i].barcode_id];
                uint8_t sl = exact_distance(b, r, UINT8_MAX);
                index[i].distance = sl;
            }

#pragma omp atomic
            d3 += std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();

            /****************************************************************************************************
             * Step 4: Sort the index by the SL-distance
             ***************************************************************************************************/

            t0 = std::chrono::high_resolution_clock::now();

            std::sort(index, index + index_size, [](const index_element& x, const index_element& y) {
                if (x.distance < y.distance)
                    return true;
                /*if (x.distance == y.distance)
                    return x.barcode_id < y.barcode_id;*/
                return false;
            });

#pragma omp atomic
            d4 += std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();

            /****************************************************************************************************
             * Select the barcode with smallest SL-distance
             ***************************************************************************************************/

            if (index[0].distance <= rejection_threshold) {
                assigned_barcode[read_id] = index[0].barcode_id;
                /*#pragma omp critical
                                {
                                    std::cout << "assign read " << read_id << " to barcode " << index[0].barcode_id << std::endl;
                                    if (read_id != index[0].barcode_id) {
                                        std::cout << "index of read " << read_id << ":" << std::endl;
                                        for(unsigned i=0; i<MAX_INDEX_SIZE; i++) {
                                            std::cout << i << ": " << index[i].barcode_id << " " << index[i].distance << std::endl;
                                        }
                                    }
                                }*/
            } else {
                //std::cout << "read " << read_id << " rejected because SL_distance to " << index[0].barcode_id << " is " << index[0].distance << std::endl;
            }
        }

        /****************************************************************************************************
         * Output debug information
         ***************************************************************************************************/

        auto t2 = std::chrono::high_resolution_clock::now();
        double total = std::chrono::duration<double, std::milli>(t2 - t0).count();
        /*std::cout << "Step 0: " << 100*d0/total << "%" << std::endl;
        std::cout << "Step 1: " << 100*d1/total << "%" << std::endl;
        std::cout << "Step 2: " << 100*d2/total << "%" << std::endl;
        std::cout << "Step 3: " << 100*d3/total << "%" << std::endl;
        std::cout << "Step 4: " << 100*d4/total << "%" << std::endl;
        std::cout << "Step 5: " << 100*d5/total << "%" << std::endl;*/

        return assigned_barcode;
    }

};

#endif //INC_2OPT_KMER_FILTERED_CALLING_HOST_V1_H
