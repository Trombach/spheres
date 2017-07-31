#include <dlib/optimization.h>
#include <vector>
#include <iostream>
#include "geometry.h"
#include "structure.h"
#include "parameter.h"
#include "stop_strategy.h"

using namespace std;

typedef dlib::matrix<double,0,1> column_vector;

class pairPotential 
{
    private:
        virtual double E (double distance) {return distance;};
        virtual double dE_dr (double distance) {return distance;};
        //double d2E_dr2 () {};

    public:
        pairPotential();
        double calcEnergy (const column_vector &v);
        const column_vector calcGradient (const column_vector &v);
        //std::vector< std::vector<double> > calcHessian ();
};

pairPotential::pairPotential (void) {}

double pairPotential::calcEnergy (const column_vector &v)
{
    structure S;
    for (long i = 0; i < v.size() / 3; ++i) 
    {
        coord3d sphere(v(3 * i), v(3 * i + 1), v(3 * i + 2));
        S.push_back(sphere);
    }
    double f(0);
    for (int i = 0; i < S.nAtoms(); ++i) {
        for (int j = i + 1; j < S.nAtoms(); ++j) {
            f += this->E (coord3d::dist (S[i],S[j]));
        }
    }
    return f;
}

const column_vector pairPotential::calcGradient (const column_vector &v)
{
    structure S;
    for (long i = 0; i < v.size() / 3; ++i) 
    {
        coord3d sphere(v(3 * i), v(3 * i + 1), v(3 * i + 2));
        S.push_back(sphere);
    }
    vector<coord3d> gradients (S.nAtoms(), coord3d());
    for (int i = 0; i < S.nAtoms(); i++)
    {
        for (int j = i + 1; j < S.nAtoms(); j++)
        {
            coord3d distanceVector = S[j]-S[i];
            double gradValue = this->dE_dr (coord3d::dist(S[i],S[j]));
            coord3d twoBodyGradient = distanceVector / distanceVector.norm() * gradValue;
            gradients[i] += twoBodyGradient;
            gradients[j] -= twoBodyGradient;
        }
    }
    column_vector df (gradients.size() * 3);
    for (vector<coord3d>::size_type i = 0; i< gradients.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            df(3 * i + j) = gradients[i][j];
        }
    }
    return df;
}









class LJ : public pairPotential 
{
    private:
        double _epsilon;
        double _rm;
        double _exp1;
        double _exp2;
        double E (double distance);
        double dE_dr (double distance);
        

    public:
        LJ() :  _epsilon(1),
                _rm(1),
                _exp1(12),
                _exp2(6)
        {}
};

double LJ::E (double distance)
{
    cout << "this is LJ" << endl;
    return (_epsilon / (_exp1/_exp2 - 1)) * ( (pow (_rm / distance, _exp1)) - (_exp1/_exp2) * (pow (_rm / distance, _exp2)) );
}

double LJ::dE_dr (double distance)
{
    cout << "this is LJ grad" << endl;
    return ( _epsilon / (_rm * (_exp1/_exp2 - 1)) ) * ( _exp1 * (pow (_rm / distance, _exp1 + 1)) - _exp1 * (pow (_rm / distance, _exp2 + 1)) );
}


structure structure::optimize_test (ostream &min, parameter<int> switches, parameter<double> opt, pairPotential *potential)
{
    min.precision(16);

    const size_t nsteps = static_cast<const size_t>(opt.get("nsteps"));
    const double stop_crit = opt.get("convergence");

    column_vector x((this->nAtoms()) * 3);
    for (int i = 0; i < this->nAtoms(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            x(3 * i + j) = (*this)[i][j];
        }
    }

    dlib::find_min(dlib::bfgs_search_strategy(),
            dlib::stop_strategy(stop_crit, nsteps),
            potential->calcEnergy(x), potential->calcGradient(x), x, -(this->nAtoms) * 1000);

}

int main ()
{
    vector<coord3d> coords;
    coords.push_back(coord3d(0,0,0));
    coords.push_back(coord3d(1,0,0));
    coords.push_back(coord3d(1,0.87,0));
    structure test(1,coords);

    column_vector x(test.nAtoms() * 3);
    for (int i = 0; i < test.nAtoms(); ++i) 
    {
        for (int j = 0; j<=2; ++j) 
        {
            x(3 * i + j) = (test)[i][j];
        }
    }

    LJ len;
    cout << len.calcEnergy(x) << endl;
    const column_vector grad = len.calcGradient(x);
    cout << grad << endl;

return 0;
}

