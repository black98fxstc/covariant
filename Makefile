CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -I. -g
LIBS = -lfftw3 -lfftw3f -lm
HEADERS = Weighty.hpp Laplacian.hpp Covariant.hpp TestData.hpp

all: laplacian covariant testdata

covariant2: main2.o
	$(CXX) $(CXXFLAGS) -o covariant2 main2.o $(LIBS)

covariant3: main3.o
	$(CXX) $(CXXFLAGS) -o covariant3 main3.o $(LIBS)

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
	rm -f main2.o main3.o covariant2 covariant3 covariant.o laplacian.o laplacian covariant testdata.o testdata

.PHONY: all clean