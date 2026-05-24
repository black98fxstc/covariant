CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -I. -g
LIBS = -lfftw3 -lfftw3f -lm
HEADERS = Dimensions.hpp Function.hpp Weighty.hpp Laplacian.hpp Covariant.hpp TestData.hpp

all: covariant laplacian testdata

testdata: testdata.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o testdata testdata.o $(LIBS)

laplacian: laplacian.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o laplacian laplacian.o $(LIBS)

covariant: covariant.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o covariant covariant.o $(LIBS)

testdata.o: $(HEADERS)
laplacian.o: $(HEADERS)
covariant.o: $(HEADERS)
main2.o: $(HEADERS)
main3.o: $(HEADERS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o laplacian covariant testdata
	rm -f *.exe
	
.PHONY: all clean