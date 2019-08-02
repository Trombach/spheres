#ifndef STORAGE
#define STORAGE

#include <map>
#include <iostream>
#include <cassert>
#include "structure.h"

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
