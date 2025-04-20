CXX       = c++
CXXFLAGS   = -O3 -Wall -shared -std=c++17 
DEVFLAGS = -Wall -shared -ggdb -O0 -std=c++17 

UNAME_S = $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    CXXFLAGS += -undefined dynamic_lookup
	DEVFLAGS += -undefined dynamic_lookup
endif

all:
	$(CXX) $(CXXFLAGS) -fPIC $(shell python3 -m pybind11 --includes)  pyValidatePath.cpp -o Grid_Path_Checker$(shell python3-config --extension-suffix)

dev:
	$(CXX) $(DEVFLAGS) -fPIC $(shell python3 -m pybind11 --includes)  pyValidatePath.cpp -o Grid_Path_Checker$(shell python3-config --extension-suffix)

