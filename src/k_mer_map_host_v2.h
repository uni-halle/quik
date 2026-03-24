//
// Created by agkya on 12.03.26.
//

#pragma once
#include <cstdint>
#include <vector>

#include "barcode.h"
#include "barcode_set.h"

namespace barcode_calling {

    template <unsigned k, unsigned PSEUDO_DISTANCE_WINDOW_SIZE>
    class k_mer_map_host_v2 {

        static constexpr unsigned k_mer_start_position_count = BARCODE_LENGTH - k + 1;
        static constexpr uint32_t k_mer_count = 1 << 2 * k; // 4^k
        static constexpr uint32_t mask = (1 << 2 * k) - 1; // bitmask with ones at its 2k least significant positions

        static_assert(2 * k < 32, "invalid value of k");
        static_assert(k < BARCODE_LENGTH);

        /****************************************************************************************
         * For each combination of k-mer x and possible start position i, we
         * maintain two lists barcode_ids[xi] and start_positions[xi] of the same length l,
         * such that for each index j between 0 and l-1, the k-mer x starts in
         * barcodes[barcode_ids[xi][j]] at start_positions[xi][j].
         ******************************************************************************************/

        std::array<std::vector<unsigned>, k_mer_count * k_mer_start_position_count> barcode_ids;
        std::array<std::vector<uint8_t>, k_mer_count * k_mer_start_position_count> start_positions;

        unsigned flat(unsigned x, unsigned i) const {
            return x * k_mer_start_position_count + i;
        }

    public:

        /**
         * Construct a k-mer-map which maintains a list map[x][i] for each combination of
         * a k_mer x in {A,C,G,T}^k and a possible start position i  i in [0...BARCODE_LENGTH-k].
         *
         * The list barcode_id[x][i] contains the ids of all barcodes in which the k-mer x
         * starts at a position j in a window of PSEUDO_DISTANCE_WINDOW_SIZE around i.
         *
         * Thus, barcodes[x][i] contains all barcodes b for which x == b[j...j+k-1] for a
         * start position j in [i-PSEUDO_DISTANCE_WINDOW_SIZE, i+PSEUDO_DISTANCE_WINDOW_SIZE].
         *
         * @param barcodes
         */
        k_mer_map_host_v2(const barcode_set& barcodes)
            : k_mer_map_host_v2(barcodes, 0, barcodes.size()){}

        /**
         * Construct a k-mer-map which maintains a list map[x][i] for each combination of
         * a k_mer x in {A,C,G,T}^k and a possible start position i  i in [0...BARCODE_LENGTH-k].
         *
         * The list barcode_id[x][i] contains the ids of all barcodes in which the k-mer x
         * starts at a position j in a window of PSEUDO_DISTANCE_WINDOW_SIZE around i.
         *
         * Thus, k_mer_map[x][i] contains all barcodes b for which x == b[j...j+k-1] for a
         * start position j in [i-PSEUDO_DISTANCE_WINDOW_SIZE, i+PSEUDO_DISTANCE_WINDOW_SIZE].
         *
         * Use only the barcodes in the given range.
         **/
        k_mer_map_host_v2(const barcode_set& barcodes,
                       const size_t barcode_start_id,
                       const size_t barcode_end_id)  {

            // for each k-mer position i in the barcode
#pragma omp parallel for
            for (unsigned i = 0; i <= BARCODE_LENGTH - k; i++) {

                // for each barcode
                for (unsigned barcode_id = barcode_start_id;
                     barcode_id < barcode_end_id && barcode_id < barcodes.size();
                     barcode_id++) {

                    const barcode& b = barcodes[barcode_id];

                    // consider the k-mer x = b[i...i+k-1]
                    uint32_t x = 0; // x = b[i...i+k-1]
                    for (unsigned l = 0; l < k; l++)
                        x = (x << 2) + b[i + l].to_uint8();
                    assert(x < k_mer_count);

                    /***********************************************************************
                     * The k-mer x starts in barcode b at position i.
                     *
                     * Thus, it occurs in all lists map[xj] for all j in the interval
                     *   [i-PSEUDO_DISTANCE_WINDOW_SIZE, i+PSEUDO_DISTANCE_WINDOW_SIZE]
                     ************************************************************************/

                    unsigned j = 0;
                    if (i > PSEUDO_DISTANCE_WINDOW_SIZE)
                        j = i - PSEUDO_DISTANCE_WINDOW_SIZE;
                    for (; j < i + PSEUDO_DISTANCE_WINDOW_SIZE && j + k <= BARCODE_LENGTH; j++) {
                        unsigned xj = flat(x, j);
                        barcode_ids[xj].push_back(barcode_id);
                        start_positions[xj].push_back(j);
                    }
                }
            }
        }

        /**
         * Return the list of barcode (ids) in which the k-mer x starts within a
         * window of PSEUDO_DISTANCE_WINDOW_SIZE around the position i.
         *
         * @param x K-mer in binary representation.
         * @param i Center of start positions.
         * @return
         */
        const std::vector<unsigned>& get_barcode_ids(unsigned x, unsigned i) const {
            unsigned xi = flat(x,i);
            return barcode_ids[xi];
        }

        /**
         * Return the k-mer start positions j associated to the barcode ids that can be retrieved
         * from get_barcode_ids.
         *
         * @param x K-mer in binary representation.
         * @param i Center of start positions.
         * @return
         */
        const std::vector<unsigned>& get_start_positions(unsigned x, unsigned i) const {
            unsigned xi = flat(x,i);
            return start_positions[xi];
        }

        /**
         * Return the number of k-mers.
         * @return
         */
        __host__ __device__ static constexpr unsigned get_k_mer_count() {
            return k_mer_count;
        }
    };
}
#endif //K_MER_MAP_V1_H
