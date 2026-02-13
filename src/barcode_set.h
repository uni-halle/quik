//
// Created by agkya on 29.01.26.
//

#ifndef BARCODE_SET_H
#define BARCODE_SET_H

#include <unordered_map>
#include <vector>
#include "barcode.h"

namespace barcode_calling {

    class barcode_set : public std::vector<barcode> {

    protected:
        // each barcode has a unique identifier
        std::vector<std::string> names;

        // to find for each name its index
        std::unordered_map<std::string, unsigned> index_of;

    public:

        void add(const barcode& b, const std::string& name) {
            this->push_back(b);
            names.push_back(name);
            index_of[name] = names.size() - 1;
        }

        const std::string& get_name_of(unsigned barcode_id) const {
            return names[barcode_id];
        }

        unsigned get_index_of(const std::string& barcode_name) const {
            if (index_of.contains(barcode_name))
                return index_of.at(barcode_name);
            return UINT_MAX;
        }

    };

}

#endif //BARCODE_SET_H
