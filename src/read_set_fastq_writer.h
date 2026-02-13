//
// Created by agkya on 29.01.26.
//

#ifndef READ_SET_FASTQ_WRITER_H
#define READ_SET_FASTQ_WRITER_H

#include "read_set.h"
#include <sstream>
#include <ostream>

namespace barcode_calling {

    class read_set_fastq_writer {

    protected:
        std::stringstream ss;

    public:

        /**
         * Load a read set from a FASTQ file.
         * @param filename
         */
        read_set_fastq_writer(const read_set& reads) {
            for (unsigned read_id = 0; read_id < reads.size(); ++read_id) {
                read r = reads[read_id];
                ss << "@" << reads.get_name_of(read_id) << '\n';
                ss << r.to_string() << '\n';
                ss << "+\n";
                for (int i = 0; i < r.length(); ++i)
                    ss << 'I';
                ss << '\n';
            }
        }

        void write(std::ostream& out) const {
            out << ss.str();
        }

    };

}

#endif //READ_SET_FASTQ_WRITER_H
