CXX = g++
FLAGS = -ggdb 

all: m+.o aStar.o 
	${CXX} ${FLAGS} -o m+ m+.o aStar.o -fopenmp

m+.o: m+.cpp m+.hpp
	${CXX} ${FLAGS} -c m+.cpp -fopenmp

aStar.o: aStar.cpp m+.hpp
	${CXX} ${FLAGS} -c aStar.cpp -fopenmp

clean:
	rm -f *.o




#original make
#all: m+.o aStar.o 
#	g++ -o m+ m+.o aStar.o -fopenmp
#
#m+.o: m+.cpp m+.hpp
#	g++ -c m+.cpp -fopenmp
#
#aStar.o: aStar.cpp m+.hpp
#	g++ -c aStar.cpp -fopenmp
#
#clean:
#	rm *.o
