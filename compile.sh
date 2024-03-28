#!/bin/bash

OUTPUT_NAME="decode_caen_coincidences"

# Compile each source file separately into object files
g++ -g -c decode_caen_coincidences.C
g++ -g -c functions/functions_visualization.C
g++ -g -c functions/functions_hexdump_processing.C
g++ -g -c functions/functions_time_analysis.C

# Link the object files together to create the final executable
g++ decode_caen_coincidences.o functions_visualization.o functions_hexdump_processing.o functions_time_analysis.o -o "$OUTPUT_NAME"

# Clean up object files
rm *.o
