CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -I. -g
LIBS = -lfftw3 -lfftw3f -lm

all: covariant2

covariant2: main2.o
	$(CXX) $(CXXFLAGS) -o covariant2 main2.o $(LIBS)

covariant3: main3.o
	$(CXX) $(CXXFLAGS) -o covariant3 main3.o $(LIBS)

main2.o: Covariant.hpp TestData.hpp
main3.o: Covariant.hpp

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f main2.o main3.o covariant2 covariant3

.PHONY: all clean