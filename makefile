CXX = g++
FLAGS = -ggdb -Wall

all: m+.o mp.o aStar.o 
	${CXX} ${FLAGS} -o m+ m+.o mp.o aStar.o -fopenmp

m+.o: m+.cpp m+.hpp
	${CXX} ${FLAGS} -c m+.cpp

mp.o: mp.cpp m+.hpp
	${CXX} ${FLAGS} -c mp.cpp -fopenmp

aStar.o: aStar.cpp m+.hpp
	${CXX} ${FLAGS} -c aStar.cpp -fopenmp

clean:
	rm -rf *.o
	




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
