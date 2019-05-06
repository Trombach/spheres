#include <iostream>
#include "storage.h"
#include "structure.h"

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
