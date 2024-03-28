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
#include<algorithm>

#include "functions_time_analysis.h"

// * Function reproducing IDQ start and stop alogorithm from two vectors containing time stamps in two 
// * channels

std::pair<std::vector<double>, std::vector<double>> correlate_idq(std::vector<double> photons_A, std::vector<double> photons_B, double binsize, int bins_count) {

    double max_delay = binsize * bins_count;
    std::vector<double> bins;

    for (int i = 0; i < bins_count; ++i) {
        bins.push_back(i * binsize);
    }
    
    int i = 0;
    int j = 0;
    int len_A = photons_A.size();
    int len_B = photons_B.size();
    int counter1 = 0;
    int counter2 = 0;
    
    std::sort(photons_A.begin(),photons_A.end());
    std::sort(photons_B.begin(),photons_B.end());

    std::vector<double> coincidences(bins_count, 0.0);
    std::vector<double> delays;


    while (true) {
        if (i == len_A - 1 or j == len_B) {
            break;
        }
        
        double t_A1 = photons_A[i];
        double t_A2 = photons_A[i + 1];
        double t_B = photons_B[j];
        
        if (t_B < t_A1) {
            j++;
            continue;
        }
        
        if (t_A1 <= t_B && t_B < t_A2) {
     
            double delay = t_B - t_A1;
            if (delay < max_delay) {
                int index = static_cast<int>(std::floor(delay / binsize));
                coincidences[index]++;
            
            }
            j++;
            counter1++;
            
        } else {
            i++;
            counter2++;
        }
    }
    return std::make_pair(bins, coincidences);;
}



std::pair<std::vector<double>, std::vector<double>> correlate_time_stamps(std::vector<double> photons_A, std::vector<double> photons_B, double binsize, int bins_count) {

    double max_delay = binsize * bins_count;
    std::vector<double> bins;

    for (int i = 0; i < bins_count; ++i) {
        bins.push_back(i * binsize);
    }
    
    std::vector<double> coincidences(bins_count, 0.0);
    
    int k = 0;

    for (int j = 0; j < photons_A.size() - 1; ++j) {
        double t_start = photons_A[j];
        double t_start_2 = photons_A[j + 1];

        while (k < photons_B.size() && photons_B[k] < t_start_2) {
        
        	  if(photons_B[k]<=t_start){
        	  k++;
        	  continue;
        	  }
        	  
            double delay = photons_B[k] - t_start;
            if (delay < max_delay) {
                int index = static_cast<int>(std::floor(delay / binsize));
                coincidences[index]++;
            }
            k++;
        } 
    }
    return std::make_pair(bins, coincidences);
}


std::pair<std::vector<double>, std::vector<double>> correlate_time_stamps_2nd_order(std::vector<double> photons_A, std::vector<double> photons_B, double binsize, int bins_count) {
    double max_delay = binsize * bins_count;

    // Calculate mean and standard deviation of photons_A
    double mean_A = 0.0, std_dev_A = 0.0;
    for (double t : photons_A) {
        mean_A += t;
    }
    mean_A /= photons_A.size();
    for (double t : photons_A) {
        std_dev_A += pow(t - mean_A, 2);
    }
    std_dev_A = sqrt(std_dev_A / photons_A.size());
    
    double mean_B = 0.0, std_dev_B = 0.0;
    for (double t : photons_B) {
        mean_B += t;
    }
    mean_B /= photons_B.size();
    for (double t : photons_B) {
        std_dev_B += pow(t - mean_B, 2);
    }
    std_dev_B = sqrt(std_dev_B / photons_B.size());


    std::vector<double> delays;
    std::vector<double> correlations;

    // Iterate over different delays
    for (int d = 0; d < bins_count; ++d) {
        double delay = d * binsize;
        double cross_corr = 0.0, sum_A = 0.0, sum_B = 0.0, sum_AB = 0.0;

        // Calculate cross-correlation at current delay
        for (double t_start : photons_A) {
            for (double t_end : photons_B) {
                if (fabs(t_end - t_start) <= delay) {
                    sum_A += (t_start - mean_A) * (t_start - mean_A);
                    sum_B += (t_end - mean_B) * (t_end - mean_B);
                    sum_AB += (t_start - mean_A) * (t_end - mean_B);
                }
            }
        }

        // Calculate Pearson correlation coefficient
        double correlation_coeff = sum_AB / (sqrt(sum_A) * sqrt(sum_B));

        // Store delay and correlation coefficient
        delays.push_back(delay);
        correlations.push_back(correlation_coeff);
    }

    return std::make_pair(delays, correlations);
}


/*std::pair<std::vector<double>, std::vector<double>> rebin_photons_stream(std::vector<double> photons_stream, double binsize){
    
    
    return std::make_pair
}*/



// * Function to compute auto-correlation of a single channel (delta t). It returns the histogram bins and heights. 

std::pair<std::vector<double>, std::vector<double>> auto_correlation(std::vector<double> photons_A, double binsize, int bins_count) {
  

    double max_delay = binsize * bins_count;
    std::vector<double> bins;

    for (int i = 0; i < bins_count; ++i) {
        bins.push_back(i * binsize);
    }
    
    int i = 0;
    
    int len_A = photons_A.size();
    int counter1 = 0;
    

    std::vector<double> coincidences(bins_count, 0.0);
   
  

    while (i<len_A-1) {

        double delta = photons_A[i+1]-photons_A[i];
      
            if (delta < max_delay) {
                int index = static_cast<int>(std::floor(delta / binsize));
                coincidences[index]++;  
            }
        i++;
       
  
     }
      

     return std::make_pair(bins, coincidences);;

}



// * Function that selcts ToAs of only a single channel from the data stream. 

std::vector<double> mask_data_channel(std::vector<double> photons_ToA, std::vector<int> channels, int target_channel){
  std::vector<double> masked_data = {};

  for (size_t k = 0; k < photons_ToA.size(); ++k) {
        if (channels[k]==target_channel) {
            masked_data.push_back(photons_ToA[k]);
        }
    }
  return masked_data;
}


