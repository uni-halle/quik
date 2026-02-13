//
// Created by agkya on 28.01.26.
//

#ifndef EXTENDED_BARCODE_ASSIGNMENT_H
#define EXTENDED_BARCODE_ASSIGNMENT_H
#include "barcode_assignment.h"

namespace barcode_calling {

    class extended_barcode_assignment : public barcode_assignment {

    protected:

        /**
         * To each read r, we assign two different barcodes:
         *
         *   1) A barcode b1 to which r has minimum distance.
         *   2) A barcode b2 to which r has second best distance.
         */
        std::vector<unsigned> barcode_id_2nd;
        std::vector<int> distance_to_2nd_barcode;

    public:

        /**
         * Define the best-fitting barcode to some read.
         * @param read_id
         * @param barcode_id
         * @param distance
         */
        void assign_as_1st_barcode(unsigned read_id, unsigned barcode_id, int distance) {
            assign_read_to_barcode(read_id, barcode_id, distance);
        }

        /**
         * Assign the second-best-fitting barcode to some read.
         * @param read_id
         * @param barcode_id
         * @param distance
         */
        void assign_as_2nd_barcode(unsigned read_id, unsigned barcode_id, int distance) {
            barcode_id_2nd[read_id] = barcode_id;
            distance_to_2nd_barcode[read_id] = distance;;
        }

    public:

        extended_barcode_assignment(unsigned read_count) :
            barcode_assignment(read_count),
            barcode_id_2nd(read_count, UINT_MAX),
            distance_to_2nd_barcode(read_count, INT_MAX) {}

        /**
        * Return an barcode index minimum distance to the read.
        * @param read_id
        * @return
        */
        unsigned get_1st_barcode(unsigned read_id) const {
            return get_assigned_barcode(read_id);
        }

        /**
         * Return an barcode index second-to-minimum distance to the read.
         * @param read_id
         * @return
         */
        unsigned get_2nd_barcode(unsigned read_id) const {
            return barcode_id_2nd[read_id];
        }

        int get_distance_to_1st_barcode(unsigned read_id) const {
            return distance[read_id];
        }

        int get_distance_to_2nd_barcode(unsigned read_id) const {
            return distance_to_2nd_barcode[read_id];
        }

    };

}

#endif //EXTENDED_BARCODE_ASSIGNMENT_H
