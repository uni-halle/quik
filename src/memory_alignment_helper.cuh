//
// Created by agkya on 25.02.26.
//

#ifndef BARCODE_CALLING_ALIGNMENT_HELPER_CUH
#define BARCODE_CALLING_ALIGNMENT_HELPER_CUH
#include <cstdint>

namespace quik {

    class memory_alignment_helper {

        unsigned char* base;
        size_t offset = 0;

        __host__ __device__
        static size_t align_up(size_t x, size_t a) {
            return (x + a - 1) & ~(a - 1);
        }

    public:

        __host__ __device__
        memory_alignment_helper(unsigned char* base) : base(base) {}

        template<typename type>
        __host__ __device__
        type* register_array(size_t count) {
            offset = align_up(offset, alignof(type));
            type *ptr = reinterpret_cast<type*>(base + offset);
            offset += count * sizeof(type);
            return ptr;
        }

        /**
         * Return the number of bytes.
         * @return
         */
        __host__ __device__
        size_t size() const {
            return offset;
        }
    };

}

#endif //BARCODE_CALLING_ALIGNMENT_HELPER_CUH