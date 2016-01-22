CC= clang++

LDFLAGS= -I/usr/local/Cellar/gsl/1.16/include -I/usr/local/Cellar/libconfig/1.5/include
LLIBFLAGS= -L/usr/local/lib
FLAGS= $(LDFLAGS) $(LLIBFLAGS)

all:
	$(CC) -O3 $(FLAGS) -lconfig++ -lgsl -lgslcblas -std=c++11 -o optimize main.cc structure-optimize.cc lina.cc structure.cc iop.cc
