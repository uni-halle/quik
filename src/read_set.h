//
// Created by agkya on 29.01.26.
//

#ifndef READ_SET_H
#define READ_SET_H

#include <unordered_map>
#include <vector>
#include "read.h"

namespace barcode_calling {

    class read_set {

    protected:

        /***********************************************************************
         * The read sequence data is concatenated into one single large array.
         *
         * Given some read_id, the associated sequence data of reads[read_id]
         * is located at the positions read_data_start_index[read_id] to
         * read_data_start_index[read_id+1] - read_data_start_index[read_id].
         ************************************************************************/

        std::vector<uint8_t> read_data;
        std::vector<size_t> read_data_start_index = {0};

        // each read has a unique identifier
        std::vector<std::string> names;

        // each read has a length (number of ASCII characters)
        std::vector<unsigned> read_sequence_lengths;

        // to find the index of a read from its name
        std::unordered_map<std::string, unsigned> index_of;

        // temporary working vector
        std::vector<uint8_t> tmp;

    public:

        read_set() = default;

        /**
         * Add another read to the set.
         * @param sequence ASCII Sequence.
         * @param fastq_id Sequence ID from the FASTQ file.
         */
        read add(const std::string& sequence, const std::string& fastq_id) {

            unsigned read_id = names.size();

            // insert a new read at the end of the list
            index_of[fastq_id] = read_id;
            names.push_back(fastq_id);
            read_sequence_lengths.push_back(sequence.length());

            /*************************************************************************
             * Convert the ASCII string into a array of bytes.
             * array_length = ceil(sequence.size()/4) as 4 bases fit into one uint8_t
             *************************************************************************/

            size_t array_length = (sequence.size() + 4 - 1) / 4;
            tmp.resize(array_length);
            sequence::string_to_bytes(sequence, tmp.data());

            // insert the converted sequence at the end of the read_data vector
            read_data_start_index.push_back(read_data_start_index.back() + array_length);
            read_data.insert(read_data.end(), tmp.begin(), tmp.end());

            return operator[](read_id);
        }

        /**
         * Return the number of reads in this set.
         * @return
         */
        size_t size() const {
            return names.size();
        }

        /**
         * Does the read set contain no reads?
         * @return
         */
        bool empty() const {
            return size() == 0;
        }

        const std::string& get_name_of(unsigned read_id) const {
            return names[read_id];
        }

        unsigned get_index_of(const std::string& read_name) const {
            if (index_of.contains(read_name))
                return index_of.at(read_name);
            return UINT_MAX;
        }

        read operator[](unsigned read_id) const {
            const uint8_t* start_ptr = read_data.data() + read_data_start_index[read_id];
            return read(start_ptr, read_sequence_lengths[read_id]);
        }

        /**
         * Return the vector of concatenated read data.
         */
        const std::vector<uint8_t>& get_sequence_data() const {
            return read_data;
        }

        const std::vector<unsigned>& get_sequence_lengths() const {
            return read_sequence_lengths;
        }

        /**
         * Return the vector of start positions.
         * @return
         */
        const std::vector<size_t>& get_start_indices() const {
            return read_data_start_index;
        }

    };

}

#endif //READ_SET_H
