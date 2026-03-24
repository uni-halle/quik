//
// Created by agkya on 29.01.26.
//

#pragma once

#include <vector>
#include "read.h"
#include "read_file.h"

namespace barcode_calling {

    class read_set {

    protected:
        /***********************************************************************
         * A read set is defined as a subset of coherent reads from a read file
         ************************************************************************/

        const size_t read_count;

        /***********************************************************************
         * The read sequence data is concatenated into one single large array.
         *
         *    read_data = uint8_t[read_data_size].
         *
         * Given some read_id, its associated sequence data are stored in the
         * array read_data beginning at position read_data_start_index[read_id].
         * Its data end at position read_data_start_index[read_id+1]. Thus,
         *
         *   read_data_start_index = size_t[read_count+1].
         *
         * Each read has a length (number of ASCII characters)
         *
         *   read_sequence_lengths = unsigned[read_count].
         *
         ************************************************************************/

        uint8_t* read_data = nullptr;
        size_t* read_data_start_index = nullptr;
        unsigned* read_lengths = nullptr;

        size_t read_data_size;

    public:
        /**
         * Create a read set of all reads included in a redd file.
         *
         * @param reads Read file containing all reads.
         */
        explicit read_set(const read_file& reads) : read_set(reads, 0, reads.size()) {}

        /**
         * Create a read set containing a coherent subset of reads from a file.
         *
         * @param reads Read file containing all reads.
         * @param read_id_begin Start id in [0,..., reads.size()].
         * @param read_id_end End id in [0,..., reads.size()] (exclusively).
         */
        read_set(const read_file& reads,
                 const size_t read_id_begin,
                 const size_t read_id_end)
            : read_count(read_id_end - read_id_begin) {
            assert(read_id_begin <= read_id_end);
            assert(read_id_end <= reads.size());

            const auto read_data_begin = reads.get_sequence_data().data()
                + reads.get_start_indices()[read_id_begin];
            const auto read_data_end = reads.get_sequence_data().data()
                + reads.get_start_indices()[read_id_end];
            read_data_size = (read_data_end - read_data_begin) * sizeof(uint8_t);

            // allocate pinned memory on host side
            CUDA_CHECK(cudaMallocHost(&read_data, read_data_size));
            CUDA_CHECK(cudaMallocHost(&read_data_start_index, (read_count+1)*sizeof(size_t)));
            CUDA_CHECK(cudaMallocHost(&read_lengths, read_count*sizeof(unsigned)));

            // copy the read data and lengths of the relevant reads to the pinned area
            memcpy(read_data, read_data_begin, read_data_size);
            memcpy(read_lengths, reads.get_sequence_lengths().data(),
                   read_count * sizeof(unsigned));

            // the indices must be corrected by the read start index
            memcpy(read_data_start_index, reads.get_start_indices().data() + read_id_begin,
                   (read_count + 1) * sizeof(size_t));
            for (size_t read_id = 0; read_id <= read_count; ++read_id)
                read_data_start_index[read_id] -= reads.get_start_indices()[read_id_begin];
        }

        /**
         * Create a read set containing a coherent subset of reads from another read set.
         *
         * @param reads Read file containing all reads.
         * @param read_id_begin Start id in [0,..., reads.size()].
         * @param read_id_end End id in [0,..., reads.size()] (exclusively).
         */
        read_set(const read_set& reads,
                 const size_t read_id_begin,
                 const size_t read_id_end)
            : read_count(read_id_end - read_id_begin) {
            assert(read_id_begin <= read_id_end);
            assert(read_id_end <= reads.size());

            const auto read_data_begin = reads.get_data()
                + reads.get_start_indices()[read_id_begin];
            const auto read_data_end = reads.get_data()
                + reads.get_start_indices()[read_id_end];
            read_data_size = (read_data_end - read_data_begin) * sizeof(uint8_t);

            // allocate pinned memory on host side
            CUDA_CHECK(cudaMallocHost(&read_data, read_data_size));
            CUDA_CHECK(cudaMallocHost(&read_data_start_index, (read_count+1)*sizeof(size_t)));
            CUDA_CHECK(cudaMallocHost(&read_lengths, read_count*sizeof(unsigned)));

            // copy the read data and lengths of the relevant reads to the pinned area
            memcpy(read_data, read_data_begin, read_data_size);
            memcpy(read_lengths, reads.get_lengths(), read_count * sizeof(unsigned));

            // the indices must be corrected by the read start index
            memcpy(read_data_start_index, reads.get_start_indices() + read_id_begin,
                   (read_count + 1) * sizeof(size_t));
            for (size_t read_id = 0; read_id <= read_count; ++read_id)
                read_data_start_index[read_id] -= reads.get_start_indices()[read_id_begin];
        }


        /**
         * Create a read set containing a subset of reads from another read set.
         *
         * @param reads Read file containing all reads.
         * @param read_ids Subset of [0,..., reads.size()-1].
         */
        read_set(const read_set& reads,
                 const std::vector<unsigned>& read_ids)
            : read_count(read_ids.size()) {

            /***************************************************************************
             * Concatentate the data for the reads in the specified subset.
             ***************************************************************************/

            std::vector<uint8_t> read_data_sub;
            std::vector<unsigned> read_lengths_sub(read_count);
            std::vector<size_t> read_start_indices_sub(read_count + 1);

            read_data_sub.reserve(reads.get_data_size());
            read_start_indices_sub[0] = 0;

            for (unsigned local_read_id = 0; local_read_id < read_count; ++local_read_id) {

                unsigned global_read_id = read_ids[local_read_id];
                read_lengths_sub[local_read_id] = reads.get_lengths()[global_read_id];

                auto global_read_index_begin = reads.get_start_indices()[global_read_id];
                auto global_read_index_end = reads.get_start_indices()[global_read_id + 1];
                auto global_read_data_size = global_read_index_end - global_read_index_begin;

                // append the read data to the end of the vector
                size_t current_size = read_data_sub.size();
                assert(current_size == read_start_indices_sub[local_read_id]);
                read_data_sub.resize(current_size + global_read_data_size);

                memcpy(read_data_sub.data() + current_size,
                       reads.get_data() + global_read_index_begin,
                       global_read_data_size);

                read_start_indices_sub[local_read_id+1] = read_data_sub.size();
            }

            read_data_size = read_data_sub.size() * sizeof(uint8_t);

            // allocate pinned memory on host side
            CUDA_CHECK(cudaMallocHost(&read_data, read_data_size));
            CUDA_CHECK(cudaMallocHost(&read_data_start_index, (read_count+1)*sizeof(size_t)));
            CUDA_CHECK(cudaMallocHost(&read_lengths, read_count*sizeof(unsigned)));

            // copy the data of the relevant reads to the pinned area
            memcpy(read_data, read_data_sub.data(), read_data_size);
            memcpy(read_lengths, read_lengths_sub.data(), read_count * sizeof(unsigned));
            memcpy(read_data_start_index, read_start_indices_sub.data(),(read_count + 1) * sizeof(size_t));
        }


        ~read_set() {
            CUDA_CHECK(cudaFreeHost(read_data));
            CUDA_CHECK(cudaFreeHost(read_data_start_index));
            CUDA_CHECK(cudaFreeHost(read_lengths));
        }

        /**
         * Return the number of reads in this set.
         * @return
         */
        [[nodiscard]]
        size_t size() const {
            return read_count;
        }

        [[nodiscard]]
        bool empty() const {
            return size() == 0;
        }

        read operator[](const size_t read_id) const {
            assert(read_id < size());
            const auto start_ptr = read_data + read_data_start_index[read_id];
            return {start_ptr, read_lengths[read_id]};
        }

        /**
         * Return a read-only view on the vector of concatenated read data.
         */
        [[nodiscard]]
        const uint8_t* get_data() const {
            return read_data;
        }

        [[nodiscard]]
        const unsigned* get_lengths() const {
            return read_lengths;
        }

        /**
         * Return a read-only view on the start positions.
         * @return
         */
        [[nodiscard]]
        const size_t* get_start_indices() const {
            return read_data_start_index;
        }

        /**
         * Return the number of bytes required to store the concatenated read data.
         * @return
         */
        [[nodiscard]]
        size_t get_data_size() const {
            return read_data_size;
        }

        /**
         * Return the number of bytes required to store the concatenated read lengths.
         * @return
         */
        [[nodiscard]]
        size_t get_lengths_size() const {
            return read_count * sizeof(unsigned);
        }

        /**
         * Return the number of bytes required to store the start indices.
         * @return
         */
        [[nodiscard]]
        size_t get_start_indices_size() const {
            return (read_count + 1) * sizeof(size_t);
        }
    };
}
