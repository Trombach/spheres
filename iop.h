#ifndef IOP
#define IOP

#include <string>
#include <vector>
#include <iostream>
#include "structure.h"
#include "parameter.h"


void xyzout (structure &outputStructure, const std::string &name = "structure.xyz");
void xyzoutall (std::vector<structure> &outputStructures, const std::string &name = "all.xyz");
void simpleout (std::vector<structure> &outputStructures, std::stringstream &ouput);

bool justempty (std::string str);
bool fexists (const std::string &fileName);
std::vector<structure> readallstruct (const std::string& fileName);
int readsettings (parameter<double> &opt, std::vector<double> &p, parameter<int> &switches, double &scalingFactor);

template <typename T> void matrixout (std::vector< std::vector<T> > &matrix, std::ostream &out = std::cout);
#endif
