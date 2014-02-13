all: m+.o aStar.o 
	g++ -o m+ m+.o aStar.o -fopenmp

m+.o: m+.cpp m+.hpp
	g++ -c m+.cpp -fopenmp

aStar.o: aStar.cpp m+.hpp
	g++ -c aStar.cpp -fopenmp

clean:
	rm *.o
