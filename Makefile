CXX=g++

LIB_FLAGS = -lblas -llapack 

CXXFLAGS = $(ARMA_INCLUDE_FLAG) $(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT) 

all: propagator

propagator: propagator.cpp
	$(CXX) $(CXXFLAGS)  -o $@  $<  $(LIB_FLAGS)


.PHONY: clean

clean:
	rm -f propagator