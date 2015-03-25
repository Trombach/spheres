all:
	clang++ -lconfig++ -lgsl -lgslcblas -std=c++11 readfile.cc spheres.cc
