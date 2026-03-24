//
// Created by agkya on 24.02.26.
//

#include "barcode.h"
#include "read_set.h"
#include "distance/weighted_levenshtein_v1.cuh"
#include "distance/weighted_sequence_levenshtein_v1.cuh"

using namespace quik;

int main(int argc, char** argv) {

    barcode b("CCGACAGTTCACTTCTGGAGCTCACATCATCCAA");
    read_set reads;
    read r = reads.add("CCGGTATCAGTTCTACTGTGCGGATGTCACATCATG");

    alignment_costs costs("/home/agkya/Schreibtisch/costs.csv");

    int32_t lev = weighted_levenshtein_v1::evaluate(b, r, costs);
    int32_t sle = weighted_sequence_levenshtein_v1::evaluate(b, r, costs);

    std::cout << "weighted_levenshtein_v1          = " << sle << std::endl;
    std::cout << "weighted_sequence_levenshtein_v1 = " << sle << std::endl;

    return 0;
}
