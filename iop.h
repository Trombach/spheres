#ifndef IOP
#define IOP

#include <string>
#include <vector>
#include "structure.h"

using namespace std;

void xyzout (structure &outputStructure, const string &name);
void simpleout (vector<structure> &outputStructures, stringstream &ouput);

bool justempty(string str);
bool fexists (const std::string &fileName);
vector<structure> readallstruct (const std::string& fileName);

#endif
