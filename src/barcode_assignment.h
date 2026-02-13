//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_BARCODE_ASSIGNMENT_H
#define BARCODE_CALLING_BARCODE_ASSIGNMENT_H

#include <vector>

namespace barcode_calling {

    struct barcode_assignment {

    protected:

        // to each read, we assign one barcode and the distance to the assigned barcode
        std::vector<unsigned> assigned_barcode;
        std::vector<int> distance;

    public:
        /**
         *  Create an empty assignment, in which each read is unassigned.
         */
        explicit barcode_assignment(unsigned read_count)
            : assigned_barcode(read_count, UINT_MAX),
              distance(read_count, INT_MAX) {};

        /**
         * Is a specific read assigned to some barcode?
         * @param read_id
         * @return
         */
        bool is_assigned_to_some_barcode(unsigned read_id) const {
            return assigned_barcode[read_id] != UINT_MAX;
        }

        unsigned get_assigned_barcode(unsigned read_id) const {
            return assigned_barcode[read_id];
        }

        int get_distance(unsigned read_id) const {
            return distance[read_id];
        }

        void assign_read_to_barcode(unsigned read_id, unsigned barcode_id, int dist) {
            assigned_barcode[read_id] = barcode_id;
            distance[read_id] = dist;
        }

    };
}

#endif //BARCODE_CALLING_BARCODE_ASSIGNMENT_H
