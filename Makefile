all:
	clang++-3.5 -lconfig++ -lgsl -lgslcblas -std=c++11 -o optimize main.cc spheres.cc
