all:
	clang++ -O3 -lconfig++ -lgsl -lgslcblas -std=c++11 -o optimize main.cc optimizer.cc hessian.cc inertia.cc xyz.cc
