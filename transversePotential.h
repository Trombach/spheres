#include "potential.h"
#include "structure.h"
#include <iostream>

class TransversePotential
{
    private:
        pairPotential* _potential;
        column_vector _vector;
        double _eval;

    public:
        TransversePotential(pairPotential &potential) :
            _vector()
        {_potential = &potential;}

        void setVector(double &eval, column_vector &v);
        void setVector(double &eval, std::vector<double> &v);
        column_vector getVector() { return _vector; }
        double get_eval() { return _eval; }

        double getEnergy(structure S) {return _potential->calcEnergy(S);}
        double getEnergy(column_vector v) {return _potential->calcEnergy(v);}
        column_vector getTrueGradient(const column_vector &v) {return _potential->calcGradient(v);}
        column_vector getTrueGradient(structure &S) {return _potential->calcGradient(S);}
        column_vector getTransverseGradient(const column_vector &v);
        column_vector getTransverseGradient(structure &S);
        
        std::vector< std::vector<double> > getTrueHessian (structure &S) { return _potential->calcHessian(S); }

        void calcTransverseDirection (structure &S);
        column_vector stepUphill (structure &S);

        structure optimize(structure &S);
};

