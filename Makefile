CC= clang++
CFLAGS= -O3 -std=c++11
LFLAGS= -lconfig++ -lgsl -lgslcblas
FLAGS= $(CFLAGS) $(LFLAGS)

all:
	$(CC) $(FLAGS) -o optimize main.cc structure-optimize.cc lina.cc structure.cc iop.cc
