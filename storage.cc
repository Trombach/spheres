#include <iostream>
#include "storage.h"
#include "structure.h"
#include "iop.h"

using namespace std;


template<>
bool Storage<compareEnergy>::addCluster (structure &S)
{
    double energy = S.getEnergy();

    typename storageMap<compareEnergy>::iterator iter = _mapping.find(energy);

    if (iter == _mapping.end())
    {
        _mapping[energy] = S;
        return true;
    }

    return false;
}

template<>
bool Storage<compareInterPartDist>::addCluster (structure &S)
{
    S.propertyInterPartDist();
    vector<double> dist = S.getInterPartDist();

    typename storageMap<compareInterPartDist>::iterator iter = _mapping.find(dist);

    if (iter == _mapping.end())
    {
        _mapping[dist] = S;
        return true;
    }

    return false;
}

template<typename T>
void Storage<T>::printKeys(ostream& out)
{
    for (auto& i : _mapping) out << i.first << endl; 
}

template<typename T>
void Storage<T>::printStructures()
{
    int n(0);
    for (auto& i : _mapping) 
    {
        xyzout (i.second, to_string(n) + ".xyz");
        n++;
    }
}

template class Storage<compareEnergy>;
template class Storage<compareInterPartDist>;
