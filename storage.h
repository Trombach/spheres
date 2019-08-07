#ifndef STORAGE
#define STORAGE

#include <map>
#include <iostream>
#include <cassert>
#include "structure.h"

/*
This class (Storage) implements a wrapper for a map container, that can be used
to store and compare optimised structures. The comparator function is the
template parameter and can be redefined by the user.

Implementations can be found in "storage.cc".

Parameters
----------

storageMap :
    map container that stores structures. Third parameter is the comparator
    functor that can be redefined. Currently, "compareEnergy" and
    "compareInterPartDist" are implemented. Comparator functors need to typedef
    a key type for the map container as "keyType".
_mapping :
    mapping container.

Methods
-------

bool addCluster (structure) :
    Checks if a cluster is already stored in the map container. If not it will
    be added. Returns true if cluster was added and false otherwise.
int getSize() :
    returns number of stored structures.
void printKeys (std::ostream) :
    print stored keys to ostream.
void printStructures() :
    uses xyzoutall function from "iop.h" to print all structures to output folder.
typename storageMap<T>::iterator begin()
typename storageMap<T>::iterator end() :
    wrapper for begin() and end() methods of map container.
*/

struct compareEnergy
{
    compareEnergy(double eps = 1e-4) : _eps(eps) {}
    using keyType = double;
    double _eps;
    bool operator() (const keyType a, const keyType b) const
    {
       return (b - a > _eps);
    }
};

struct compareInterPartDist
{
    compareInterPartDist (double eps = 1e-4) : _eps(eps) {}
    using keyType = std::vector<double>;
    double _eps;
    bool operator() (const keyType a, const keyType b) const
    {
        assert(a.size() == b.size());
        for (int i = 0; i < a.size(); i++)
        {
            if (fabs(b[i] - a[i]) < _eps) continue;
            if (b[i] - a[i] > _eps) return true;
            if (b[i] - a[i] < - _eps) return false;
        }
        return false;
    }
}; 


template <typename T> 
class Storage
{
    private:
        template <typename Q> using storageMap = std::map<typename Q::keyType, structure, Q>;
        storageMap<T> _mapping;

    public:
        Storage() : _mapping() {}
        bool addCluster (structure &S);
        int getSize() {return _mapping.size();}
        void printKeys(std::ostream& out);
        void printStructures();

        typename storageMap<T>::iterator begin() {return _mapping.begin();}
        typename storageMap<T>::iterator end() {return _mapping.end();}
};

using StorageByEnergy = Storage<compareEnergy>;
using StorageByInterPartDist = Storage<compareInterPartDist>;

#endif
