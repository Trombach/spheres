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
        virtual double E (double distance) = 0;
        virtual double dE_dr (double distance) = 0;
        virtual double d2E_dr2 (double distance) = 0;

    public:
        virtual ~pairPotential() {};
        double calcEnergy (const column_vector &v);
        double calcEnergy (structure &S);
        const column_vector calcGradient (const column_vector &v);
        const column_vector calcGradient (structure &S);
        std::vector< std::vector<double> > calcHessian (structure &S);
        
        structure optimize (std::ostream &min, structure &S, parameter<int> &switches, parameter<double> &opt);

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

        LJ (double epsilon, double rm, double exp1, double exp2) : _epsilon(epsilon),
                                                                   _rm(rm),
                                                                   _exp1(exp1),
                                                                   _exp2(exp2)
        {}

        static LJ *readPotential ();
};

/*--------------------------------------------------------------------------------------*/
//                          extended LJ potential derived class
/*--------------------------------------------------------------------------------------*/

class ELJ : public pairPotential
{
    private:
        std::vector<double> _c;
        double E (double distance);
        double dE_dr (double distance);
        double d2E_dr2 (double distance);

    public:
        ELJ() : _c(30)
        {}

        ELJ(std::vector<double> c) : _c(c)
        {} 

        static ELJ *readPotential ();
};


/*--------------------------------------------------------------------------------------*/
//                          LJ potential with range cutoff derived class
/*--------------------------------------------------------------------------------------*/

class RangeLJ : public pairPotential
{
    private:
        double _epsilon;
        double _rm;
        double _exp1;
        double _exp2;
        double _range; //in units of _rm
        double E (double distance);
        double dE_dr (double distance);
        double d2E_dr2 (double distance);

    public:
        RangeLJ() : _epsilon(1),
                    _rm(1),
                    _exp1(12),
                    _exp2(6),
                    _range(1)
        {}

        RangeLJ (double epsilon, double rm, double exp1, double exp2, double range) :   _epsilon(epsilon),
                                                                                        _rm(rm),
                                                                                        _exp1(exp1),
                                                                                        _exp2(exp2),
                                                                                        _range(range)
        {}

        static RangeLJ *readPotential ();
};
#endif
