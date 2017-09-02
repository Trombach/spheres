#ifndef POTENTIAL
#define POTENTIAL

#include <dlib/optimization.h>
#include "structure.h"


typedef dlib::matrix<double,0,1> column_vector;

/*--------------------------------------------------------------------------------------*/
//                          pair potential base class
/*--------------------------------------------------------------------------------------*/

class pairPotential 
{
    private:
        virtual double E (double distance) {return distance;};
        virtual double dE_dr (double distance) {return distance;};
        virtual double d2E_dr2 (double distance) {return distance;};

    public:
        pairPotential();
        double calcEnergy (const column_vector &v);
        const column_vector calcGradient (const column_vector &v);
        std::vector< std::vector<double> > calcHessian (structure &S);
        
        template <class T> structure optimize (std::ostream &min, structure &S, parameter<int> &switches, parameter<double> &opt);
};


/*--------------------------------------------------------------------------------------*/
//                          LJ potential derived class
/*--------------------------------------------------------------------------------------*/

class LJ : public pairPotential 
{
    private:
        double _epsilon;
        double _rm;
        double _exp1;
        double _exp2;
        double E (double distance);
        double dE_dr (double distance);
        double d2E_dr2 (double distance);

    public:
        LJ() :  _epsilon(1),
                _rm(1),
                _exp1(12),
                _exp2(6)
        {}

        LJ (double epsilon = 1, double rm = 1, double exp1 = 12, double exp2 = 6) : _epsilon(epsilon),
                                                                                    _rm(rm),
                                                                                    _exp1(exp1),
                                                                                    _exp2(exp2)
        {}
};

#endif
