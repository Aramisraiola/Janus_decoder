#ifndef FUNCTIONS_TIME_ANALYSIS_H_INCLUDED
#define FUNCTIONS_TIME_ANALYSIS_H_INCLUDED


std::pair<std::vector<double>, std::vector<double>> correlate_idq(std::vector<double> photons_A, std::vector<double> photons_B, double binsize, int bins_count); 
std::pair<std::vector<double>, std::vector<double>> auto_correlation(std::vector<double> photons_A, double binsize, int bins_count);
std::vector<double> mask_data_channel(std::vector<double> photons_ToA, std::vector<int> channels, int target_channel);
std::pair<std::vector<double>, std::vector<double>> correlate_time_stamps(std::vector<double> photons_A, std::vector<double> photons_B, double binsize, int bins_count);
std::pair<std::vector<double>, std::vector<double>> correlate_time_stamps_2nd_order(std::vector<double> photons_A, std::vector<double> photons_B, double binsize, int bins_count);


#endif
