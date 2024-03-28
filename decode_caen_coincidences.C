/********************************************************************************

* FILE_NAME: decode_caen_coincidences.cpp

* AUTHOR: QUASAR project

* DESCRIPTION: Decoder for caen FERS-5200 TDC output file (.dat) from Janus. 
  Takes RunXX_list.dat file and returns a .txt file containing the histogram
  (time bin and height information) of the coincidences between two channels. 
  
* NOTES: The CAEN TDC output is organized in events that are sent out each fixed
  amount of time. Each event is a packets of the hits (detections) that took
  place in the event's time window. I kept the same notation that was used in
  the documentation for hits, events and Time of Arrival (ToA), which represents
  the hit time. 
  
********************************************************************************/


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


#include "functions/functions_hexdump_processing.h"
#include "functions/functions_time_analysis.h"
#include "functions/functions_visualization.h"

// Decoding parameters (DO NOT CHANGE)
#define _HEADER_SIZE 33
#define _EVENT_SIZE_LENGTH 2
#define _BITS_PER_HEX 8
#define _HIT_ToA_SIZE 11
const double LSB = 3.125;



int main(int argc,char **argv){

  //  ------- 1 FILE OPENING AND INITIALIZATION ------- 

  // * Declare and open file to read (pointer)
  FILE* ptr =NULL;
  ptr = fopen(argv[1],"rb");

  // * Check if file opened correctly
  if (ptr==NULL){
    std::cout<<"It was not possible to open the file"<<argv[1]<<std::endl;
    return 1;
  }

  if (argc!=4){
    std::cout<<"3 inputs must be passed in command line: .dat filename, channel start and channel stop."<<std::endl;
    return 0;
  }
  
  // * Declaration of first array of chars containing hex from header (TOT: 33 bytes)
  unsigned char header[_HEADER_SIZE];
  int r;

  int channel_start = std::atoi(argv[2]);
  int channel_stop = std::atoi(argv[3]);

  int bins_tot = 20000;

  // * Set file position indicator to 0 to start the stream
  std::fseek(ptr, 0, SEEK_END);

  // * Determination of file expected size in bit
  long file_size = std::ftell(ptr);

  // * Print file size
  std::cout<<"-------------  OPENING FILE  -------------"<<"\n"<<std::endl;
  std::cout<<"File size: " << std::setprecision(5) << (double) file_size*1E-9<<" Gb"<<"\n"<<std::endl;

  // * Close file and reopen it to reset stream to 0
  fclose(ptr);
  ptr = fopen(argv[1],"rb");

  // * Extracting header and checking again the opening
  r=fread(&header,sizeof(header),1,ptr);
  
  if (r != 1) {
      std::cout << "An error occurred while reading the file" << std::endl;
      fclose(ptr);
      return 1;
    }

  // * Declaration of needed variables 
  long long event_size_sum = _HEADER_SIZE;
  unsigned int event_size[_EVENT_SIZE_LENGTH];

  int i=0;

  unsigned int board;
  unsigned int channel;
  unsigned int offset;
  
  double temporary_ToA;
  double ToA;
  double trigger_time_stamp;

  unsigned int acquisition_mode;
  unsigned int measurement_mode;
  unsigned int time_unit;

  unsigned int initial_time;

  double hits_number_START=0;
  double hits_number_STOP=0;

  bool reminder_ToA_bool=false;
  double reminder_ToA;
 

  // * Sanity check (Acquisition mode)
  acquisition_mode=(header[10] << 8) | header[9];
  measurement_mode=header[11];
  time_unit=header[12];

  sanity_check(acquisition_mode, measurement_mode, time_unit);

  std::cout<<"Start channel: " << channel_start<<std::endl;
  std::cout<<"Stop channel: " << channel_stop<<"\n"<<std::endl;

  std::cout<<"------------------------------------------"<<"\n"<<std::endl;
  std::cout<<"Unpacking file:"<<std::endl;

  // ------- 2. START OF WHILE LOOP: EVENT BY EVENT OPENING (until end of file: eof) ------- 

  std::vector<double> counts_tot(bins_tot,0);
  std::vector<double> time_bins(bins_tot,0);

  std::vector<double> photons_start_ToA = {};
  std::vector<double> photons_stop_ToA = {};

  std::vector<double> photons_start_ToA_clustered = {};
  std::vector<double> photons_stop_ToA_clustered = {};

  int events_cluster_counter=0;
  int clustering_factor = 100;
  
  while(!feof(ptr)){
    
    // * Read event size (2 bytes) and update number of bytes read in event_size_sum (+2 bytes)
    r=fread(&event_size, _EVENT_SIZE_LENGTH, 1, ptr);
    event_size_sum+=_EVENT_SIZE_LENGTH;

    // * Determine event size by swapping endianess
    unsigned int event_size_bytes = (event_size[1] << 8) | event_size[0];
    unsigned char my_event[event_size_bytes-_EVENT_SIZE_LENGTH];

    // * Declaration of 2 buffers to store event's hits information (ToA and channel) (this makes the loop 3.5 times faster!)
    const size_t buffer_size = event_size_bytes;
    int expected_size = (sizeof(my_event) - 11) / _HIT_ToA_SIZE;

    std::vector<double> buffer_ToA;
    std::vector<int> buffer_channel;
    
    buffer_ToA.resize(expected_size);;
    buffer_channel.resize(expected_size);
    

    // * Read whole event starting from byte following event size (#hits * 8 bytes)
    r=fread(&my_event, event_size_bytes-_EVENT_SIZE_LENGTH,1,ptr);
    

    // * Read event time stamp from my_event array containing whole event (without 2 first bytes from event size)
    trigger_time_stamp = U64_to_double_converter(combine_hex_little_endian(extract_subarray(my_event,0,7), _BITS_PER_HEX));
    
    // * Check if the number of bits of the event that charachterize the hits is divisible by 11 
    // * (1 byte (Board) + 1byte (Channel) + 1 byte (Edge) + 8 bytes (ToA))
    if ((sizeof(my_event)-10)%_HIT_ToA_SIZE != 0){
      std::cout<<"Problem with data unpacking, non uniform size of events";
      exit(1);
    }

    // * Display progress bar (*** --> TO DEBUG ***)
    progress_bar(event_size_sum,file_size);

    // ------- 2. START OF FOR LOOP: HIT BY HIT OPENING (UNPACKING SINGLE EVENT) -------
    
    for (int j=0; j<(sizeof(my_event)-11)/_HIT_ToA_SIZE;j++){
    
      // * Distance (offset) between same information (Board, Channel, Edge, ToA) in different hits
      offset=11*j;

      // * Read hit charachterizing parameters from my_event array 
      board = my_event[10+offset];
      channel = my_event[11+offset];

      // * Read temporary ToA in big endian and then transform it into ns depending if it was registered as LSB or as ns
      unsigned char *temporary_ToA_big_endian=extract_subarray(my_event,13+offset,20+offset);

      
      // * If statement to chose the evaluation method for the ToA, depending on its unit (LSB or ns)
      if (time_unit==0x00){
      temporary_ToA=combine_hex_little_endian(temporary_ToA_big_endian, _BITS_PER_HEX);
      ToA = temporary_ToA*LSB;
      
      }

      else if (time_unit==0x01){
        temporary_ToA=U64_to_double_converter(combine_hex_little_endian(temporary_ToA_big_endian, _BITS_PER_HEX));
        ToA = temporary_ToA*1E3;
        
        
      }

      if (i==0 && j==0){
        initial_time=ToA;
      }
      
      // * Store ToA information in buffer to write on file
      buffer_ToA[j] = ToA;  
      buffer_channel[j] = channel;
      
      // * Erase memory allocated for hit ToA 
      delete[] temporary_ToA_big_endian; 
      }
      
      
	    // * Declare photon vector corresponding to each channel and filling them with the masking function
      
      
      photons_start_ToA = mask_data_channel(buffer_ToA, buffer_channel, channel_start);
      photons_stop_ToA = mask_data_channel(buffer_ToA, buffer_channel, channel_stop);

      
      
      
      if (events_cluster_counter!=0){
        
        photons_start_ToA_clustered.insert( photons_start_ToA_clustered.end(), photons_start_ToA.begin(), photons_start_ToA.end() );
        photons_stop_ToA_clustered.insert( photons_stop_ToA_clustered.end(), photons_stop_ToA.begin(), photons_stop_ToA.end() );
  

      }
      
      else{
        photons_start_ToA_clustered = photons_start_ToA;
        photons_stop_ToA_clustered = photons_stop_ToA;

      }

      events_cluster_counter=events_cluster_counter+1;

      if (events_cluster_counter==clustering_factor){

        events_cluster_counter=0;

        std::pair<std::vector<double>, std::vector<double>> pair;
        std::vector<double> counts;

        if (photons_start_ToA_clustered.size()<=0 || photons_stop_ToA_clustered.size()<=0){
        i++;
        event_size_sum+=event_size_bytes;
        continue;
        }

        else{
          if(channel_start!=channel_stop){
            pair= correlate_idq(photons_start_ToA_clustered,photons_stop_ToA_clustered,3.125,bins_tot);
            counts = pair.second;
          }
          else{
            pair= auto_correlation(photons_start_ToA_clustered,3.125,bins_tot);
            counts = pair.second;
          }
         
        
        


      // HERE I LOOSE SOMETHING
        if (i<=clustering_factor){
        time_bins = pair.first;
        counts_tot=counts;
        
        }
      
        else {
          counts_tot = addVectors(counts_tot, counts);
        }
      

        photons_start_ToA_clustered.clear();
        photons_start_ToA_clustered.shrink_to_fit();

        photons_stop_ToA_clustered.clear();
        photons_stop_ToA_clustered.shrink_to_fit();
        }
      }
  

    // Update while loop index i and number of bytes read
    event_size_sum+=event_size_bytes; 
    i++;

    // CLose all vectors to save memory
    buffer_channel.clear(),
    buffer_channel.shrink_to_fit();

    // * Compute number of hits
    hits_number_START=hits_number_START+photons_start_ToA.size();
    hits_number_STOP=hits_number_STOP+photons_stop_ToA.size();


    buffer_ToA.clear();
    buffer_ToA.shrink_to_fit();

    photons_start_ToA.clear();
    photons_start_ToA.shrink_to_fit();

    photons_stop_ToA.clear();
    photons_stop_ToA.shrink_to_fit();
    
    //memset(my_event, 0, sizeof(my_event));
    }
    
  
  // * Print number of bytes read
  
 std::cout<<"Number of bytes decoded: "<< std::setprecision(5) << (double) event_size_sum*1E-9<<" Gb"<<"\n"<<std::endl;

 double measurement_time = (ToA-initial_time)*1E-12;
 double measurement_time_min = measurement_time / 60;
 
 std::cout<<"Time of acquisition decoded: "<< std::setprecision(3) << measurement_time_min<<" minutes"<<"\n"<<std::endl;

 std::cout<<"Number of START hits registered: "<< hits_number_START <<std::endl;
 std::cout<<"Number of STOP hits registered: "<< hits_number_STOP <<"\n"<<std::endl;

  // * Creation of output file name
  std::string filename = "output_START_" + std::to_string(channel_start) + "_STOP_" + std::to_string(channel_stop) + "_TIME_" + std::to_string(measurement_time_min) +"min.txt";

  // * Open output file and fill it with the time bins and latest value for hisogram heights
  std::ofstream outputFile(filename);

 

  if (outputFile.is_open()) {
    outputFile << "Time" << "\t" << "Coincidences"<<"\n";
        
    for (size_t i = 0; i < time_bins.size(); ++i) {
            outputFile << time_bins[i] << "\t" << counts_tot[i]<<"\n"; //<<"\t" << counts_tot_ac_A[i]<<"\t"<< counts_tot_ac_B[i]<<"\n";
        }
        outputFile.close();
    } 
    
    else {
        std::cout << "Unable to create the file." << std::endl;
    }

	
  // * File creation confirmation
  std::cout << "Histogram file created: " + filename << std::endl;
 
  return 0;
}
