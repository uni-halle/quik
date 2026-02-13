//
// Created by agkya on 27.01.26.
//

#ifndef DETERMINE_BARCODE_LENGTH_H
#define DETERMINE_BARCODE_LENGTH_H
#include <fstream>
#include <string>

/*

inline unsigned determine_barcode_length(const std::string& barcode_file) {

    // read one line of the barcode file to find the barcode length
    std::ifstream infile(barcode_file);
    std::string line;
    if (!std::getline(infile, line))
        throw std::runtime_error("Error parsing barcode file");
    infile.close();

    // extract the barcode length (number of characters)
    const unsigned BARCODE_LENGTH = line.length();

    return BARCODE_LENGTH;

}

*/

#endif //DETERMINE_BARCODE_LENGTH_H
