// Libraries import 
#include "assert.h"
#include <iostream>
#include <string>
#include <bitset>
#include <stdio.h>
#include <climits>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <utility>

#include "functions_hexdump_processing.h"

// * Function to combine a char array of hex (8-bits each) into an uint64_t with 
// * little endian endianess.

uint64_t combine_hex_little_endian(const unsigned char* hexValues, size_t size) {
    uint64_t combined = 0;

    for (size_t l = size; l > 0; --l) {
         combined = (combined << 8) | hexValues[l - 1];
    }

    return combined;
}


// * Function transforming an uint64_t into a double

double U64_to_double_converter(uint64_t val){
  double convertedValue = 0.0;
  memcpy(&convertedValue, &val, sizeof(convertedValue));
  return convertedValue;
}


// * Function transforming a double into an uint64_t

uint64_t double_to_U64_converter(double val){
  uint64_t convertedValue = 0.0;
  memcpy(&convertedValue, &val, sizeof(convertedValue));
  return convertedValue;
}


// * Function that returns a char array composed by a fraction of an original
// * array, with the elements from start to stop. Equivalent of array[start:stop]
// * in python

unsigned char* extract_subarray(unsigned char* arr, int start, int end) {
    int subarray_size = end - start + 1;
    unsigned char* subarray = new unsigned char[subarray_size];

    for (int k = 0; k < subarray_size; ++k) {
        subarray[k] = arr[start + k];
    }

    unsigned char* temp = new unsigned char[subarray_size];

    std::copy(subarray, subarray + subarray_size, temp);

    delete[] subarray; 


    return temp;
}

std::vector<double> addVectors(const std::vector<double>& vec1, const std::vector<double>& vec2) {
   
    if (vec1.size() != vec2.size()) {
        std::cout << "Vectors must be of equal size." << std::endl;
        return std::vector<double>(); }

    std::vector<double> result;
    result.reserve(vec1.size()); 

    
    for (size_t i = 0; i < vec1.size(); ++i) {
        result.push_back(vec1[i] + vec2[i]);
    }

    return result;
}



