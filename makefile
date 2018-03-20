CXX = g++-4.9
FLAGS = -std=gnu++11 -O2 -Wall

all: m+.o mp.o aStar.o 
	${CXX} ${FLAGS} -o m+1 m+.o mp.o aStar.o -fopenmp

m+.o: m+.cpp m+.hpp
	${CXX} ${FLAGS} -c m+.cpp

mp.o: mp.cpp m+.hpp
	${CXX} ${FLAGS} -c mp.cpp -fopenmp

aStar.o: aStar.cpp m+.hpp
	${CXX} ${FLAGS} -c aStar.cpp -fopenmp

clean:
	rm -rf *.o
