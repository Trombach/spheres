all:
	clang++ -lconfig++ -lgsl -lgslcblas -std=c++11 main.cc spheres.cc
