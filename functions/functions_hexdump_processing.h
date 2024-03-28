#ifndef FUNCTIONS_HEXDUMP_PROCESSING_H_INCLUDED
#define FUNCTIONS_HEXDUMP_PROCESSING_H_INCLUDED


uint64_t combine_hex_little_endian(const unsigned char* hexValues, size_t size);
double U64_to_double_converter(uint64_t val);
uint64_t double_to_U64_converter(double val);
unsigned char* extract_subarray(unsigned char* arr, int start, int end);
std::vector<double> addVectors(const std::vector<double>& vec1, const std::vector<double>& vec2);


#endif
