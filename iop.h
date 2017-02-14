#ifndef IOP
#define IOP

#include <string>
#include <vector>
#include "structure.h"
#include "parameter.h"


void xyzout (structure &outputStructure, const std::string &name);
void simpleout (std::vector<structure> &outputStructures, std::stringstream &ouput);

bool justempty(std::string str);
bool fexists (const std::string &fileName);
std::vector<structure> readallstruct (const std::string& fileName);
int readsettings (parameter<double> &opt, std::vector<double> &p, parameter<int> &switches, double &scalingFactor);
#endif
