//
// Created by agkya on 11.07.25.
//

#ifndef K_MER_MAP_V1_H
#define K_MER_MAP_V1_H
#include <cstdint>
#include <vector>

#include "../barcode.h"

namespace barcode_calling {

    template <unsigned k>
    class k_mer_map_host_v1 : public std::vector<std::vector<std::vector<unsigned>>> {

        static constexpr uint32_t map_size = 1 << 2 * k; // 4^k
        static constexpr uint32_t mask = (1 << 2 * k) - 1; // bitmask with ones at its 2k least significant positions

    public:
        /***************************************************************************************************
         * Construct a k-mer-map which contains a list of barcode ids for each combination of a
         * k-mer x in {A,C,G,T}^k and each possible k-mer start position i in [0...BARCODE_LENGTH-k].
         *
         * Thus, k_mer_map[x][i] is a list of exactly those barcode id's b, for which x == b[i...i+k-1].
         ***************************************************************************************************/

        k_mer_map_host_v1(const barcode_set& barcodes) {

            static_assert(2 * k < 32, "invalid value of k");
            static_assert(k < BARCODE_LENGTH);

            const unsigned k_mer_start_position_count = BARCODE_LENGTH - k + 1;

            // create empty k-mer-map
            resize(map_size, std::vector<std::vector<unsigned>>(k_mer_start_position_count));

            // for each k-mer position i in the barcode
#pragma omp parallel for
            for (unsigned i = 0; i <= BARCODE_LENGTH - k; i++) {

                // for each barcode
                for (unsigned barcode_id = 0; barcode_id < barcodes.size(); barcode_id++) {
                    const barcode& b = barcodes[barcode_id];

                    // consider the k-mer x = b[i...i+k-1]
                    uint32_t x = 0; // x = b[i...i+k-1]
                    for (unsigned l = 0; l < k; l++)
                        x = (x << 2) + b[i + l].to_uint8();
                    assert(x < map_size);

                    // k-mer x is located in barcode b at position i
                    (*this)[x][i].push_back(barcode_id);
                }
            }
        }

        /**
         * Return the number of k-mers.
         * @return
         */
        __host__ __device__ static unsigned get_k_mer_count() {
            return map_size;
        }
    };
}
#endif //K_MER_MAP_V1_H
