//
// Created by steffen on 24.07.24.
//

#ifndef BARCODE_CALLING_SCRAMBLE_H
#define BARCODE_CALLING_SCRAMBLE_H

#include "read.h"
#include "barcode.h"
#include <random>

std::random_device rd();

//std::default_random_engine eng(rd());
std::default_random_engine eng(0);
std::uniform_real_distribution<double> distr(0, 1);

// simulate a synthesis error
read scramble(const barcode &b, double p_ins, double p_del, double p_sub) {

    assert(p_ins >= 0 && p_ins <= 1);
    assert(p_del >= 0 && p_del <= 1);
    assert(p_sub >= 0 && p_sub <= 1);
    assert(p_ins + p_del + p_sub <= 1);

    std::stringstream ss;
    unsigned len = 0;

    char nucleotides[] = {'A', 'C', 'G', 'T'};

    // for each position in the barcode
    for (unsigned i = 0; i < barcode::LENGTH; i++) {

        double p = distr(eng);

        // which type of error occurs?
        if (p < p_ins) { // insertion

            // insert a random base
            ss << nucleotides[(unsigned) (distr(eng) * 4)];
            len++;

            i--; // b[i] is considered again in the next iteration

        } else if (p < p_ins + p_del) { // deletion
            // skip b[i]
        } else if (p < p_ins + p_del + p_sub) { // substitution

            // sample a random base != b[i]
            unsigned base;
            do {
                base = (unsigned) (distr(eng) * 4);
            } while (base == b[i]);

            ss << nucleotides[base];
            len++;
        } else { // no error
            ss << nucleotides[b[i]];
            len++;
        }
    }

    // insert random bases until the read length matches the barcode length
    for (; len < barcode::LENGTH; len++)
        ss << nucleotides[(unsigned) (distr(eng) * 4)];

    std::string str = ss.str();
    assert(str.length() <= len);

    return str.substr(0, sequence::LENGTH); // truncate additional characters
}

#endif //BARCODE_CALLING_SCRAMBLE_H
