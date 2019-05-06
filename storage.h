#ifndef STORAGE
#define STORAGE

#include <map>
#include "structure.h"

struct compareDouble
{
    compareDouble(double eps = 1e-2) : _eps(eps) {}
    double _eps;
    bool operator() (const double a, const double b) const
    {
       if (b - a > _eps) return true; 
       return false;
    }
};

class Storage
{
    private:
        typedef std::map<double, structure, compareDouble> energyStructureMap;
        energyStructureMap _mapping;

    public:
        Storage() : _mapping() {}
        bool addCluster (structure &S);
        int getSize() {return _mapping.size();}

        energyStructureMap::iterator begin() {return _mapping.begin();}
        energyStructureMap::iterator end() {return _mapping.end();}
};
#endif
