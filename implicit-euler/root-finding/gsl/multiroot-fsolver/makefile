GSL_LIB_DIR = $(GSL_HOME)/lib
GSL_INC_DIR = $(GSL_HOME)/include

OPT_FLAGS = -O3 -Wall

all: main

main: main.cc
	$(CXX) -o $@ $< $(OPT_FLAGS) -I$(GSL_INC_DIR) -L$(GSL_LIB_DIR) -lgsl -lgslcblas

clean:
	$(RM) -f main
