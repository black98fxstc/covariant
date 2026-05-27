CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -I. -g
LIBS = -lfftw3f_threads -lfftw3f -lfftw3 -lm
HEADERS = Dimensions.hpp Function.hpp Weighty.hpp Laplace.hpp Covariant.hpp TestData.hpp

all: riemann laplace testdata

testdata: testdata.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o testdata testdata.o $(LIBS)

laplace: laplace.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o laplace laplace.o $(LIBS)

riemann: riemann.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o riemann riemann.o $(LIBS)

testdata.o: $(HEADERS)
laplace.o: $(HEADERS)
riemann.o: $(HEADERS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o laplace riemann testdata
	rm -f *.exe
	
.PHONY: all clean