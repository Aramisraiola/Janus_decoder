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

#include "functions_visualization.h"

// * Function that prints a progress bar to visualize loop progression

void progress_bar(long long iteration, long total_iterations) {
    const int barWidth = 50;
    float percent = (float)iteration / total_iterations;
    int progressLength = barWidth * percent;

    std::cout << "[";
    for (int m = 0; m < barWidth; ++m) {
        if (m < progressLength)
            std::cout << "=";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(percent * 100.0) << "%\r";
    std::cout.flush();
}
void print_progress_bar(double percentage) {

    #define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
    #define PBWIDTH 60

    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
    


// * Function to check that the file being analysed was taken in the standard modes (STREAMING, LEAD_ONLY)
// * It also verifies the time unit (LSB or ns)

void sanity_check(unsigned int acquisition, unsigned int measurement, unsigned int time_stamp_unit){

  if(acquisition==0x22){
    std::cout<<"Acquisition mode: STREAMING"<<"\n"<<std::endl;
  }
  else{
    std::cout<<"Acquisition Mode is not STREAMING"<<"\n"<<std::endl;
  }
  
  if(measurement==0x01){
    std::cout<<"Measurement Mode: LEAD_ONLY"<<"\n"<<std::endl;
  }
  else{
    std::cout<<"Measurement Mode is not LEAD_ONLY"<<"\n"<<std::endl;
  }

  if (time_stamp_unit==0x00){
    std::cout<<"Time stamps are given in LSB (3.125 ps)"<<"\n"<<std::endl;
  }
  else if (time_stamp_unit==0x01){
    std::cout<<"Time stamps are given in ns"<<"\n"<<std::endl;
  }
  else{
    std::cout<<"Unrecognized time unit"<<"\n"<<std::endl;
  }
}



