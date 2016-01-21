all:
	clang++ -O3 -lconfig++ -lgsl -lgslcblas -std=c++11 -o optimize main.cc structure-optimize.cc lina.cc structure.cc iop.cc
