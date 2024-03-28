# Define compiler for whole compilation
CXX := g++

# Define additional flags
CXXFLAGS := -g -O3

# Define output program name (executable)
OUTPUT_NAME := decode_caen_coincidences

# List of source file to link
SRCS := decode_caen_coincidences.C \
        functions/functions_visualization.C \
        functions/functions_hexdump_processing.C \
        functions/functions_time_analysis.C

# Object files
OBJS := $(SRCS:.C=.o)

# Targets
all: $(OUTPUT_NAME)

$(OUTPUT_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(OUTPUT_NAME)

%.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(OUTPUT_NAME)
