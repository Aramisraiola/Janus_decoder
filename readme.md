# CAEN decoder for HBT interferometry
<div style="text-align: justify;">
This decoder is designed to interpret Janus 5203 binary output files generated by CAEN TDC (A5203) when operated in streaming mode. It specifically reads the trigger times of hits, also known as Time of Arrivals (ToA), from Janus' output files (.dat). Subsequently, it applies a coincidence algorithm to compare the start and stop times of two channels, allowing for coincidence measurements. The output of the program is a .txt file that comprises the final histogram illustrating coincidences versus time lag. The histogram is automatically binned with the finest time resolution possible, in nanoseconds (3.125 ps).
</div>

### Running the program

1. **Compiling the program:**
    The program must be compiled through the make file in order to correctly link the source codes (but can also be done through the compile.sh file)
    ```bash
     make
    ```

2. **Run the program by passing the binary caen_data.dat file, the start channel and the stop channel:**

    ```bash
   	./decode_caen_coincidences /path/to/caen_data.dat ch_start ch_stop
    ```
### Notes

- The output file is named as "output_START_{ch_start}_STOP_{ch_stop}_TIME_{acquisition time}min.txt";
- In case ch_start==ch_stop, then the autocorrelation of the target channel will be computed

