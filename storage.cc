#include <iostream>
#include "storage.h"
#include "structure.h"
#include "iop.h"

using namespace std;

bool Storage::addCluster (structure &S)
{
    double energy = S.getEnergy();

    energyStructureMap::iterator iter = _mapping.find(energy);

    if (iter == _mapping.end())
    {
        _mapping[energy] = S;
        return true;
    }

    return false;
}

void Storage::printEnergies(ostream& out)
{
    for (auto& i : _mapping) out << i.first << endl; 
}

void Storage::printStructures()
{
    int n(0);
    for (auto& i : _mapping) 
    {
        xyzout (i.second, to_string(n) + ".xyz");
        n++;
    }
}
