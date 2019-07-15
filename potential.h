#ifndef POTENTIAL
#define POTENTIAL

#include <dlib/optimization.h>
#include "stop_strategy.h"
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

    protected:
        pairPotential (
                        const int algo_switch = 1, 
                        const double stop_crit = 1e-15,
                        const int nsteps = 1000         ) 
                        : 
                        _algo_switch(algo_switch), 
                        _stop_crit(stop_crit),
                        _nsteps(nsteps) {}
        const int _algo_switch;
        const double _stop_crit;
        const int _nsteps;

    public:
        virtual ~pairPotential() {};
        double calcEnergy (const column_vector &v);
        double calcEnergy (structure &S);
        const column_vector calcGradient (const column_vector &v);
        const column_vector calcGradient (structure &S);
        std::vector< std::vector<double> > calcHessian (structure &S);

        const int getAlgoSwitch() const {return _algo_switch;}
        const double getStopCrit() const {return _stop_crit;}
        const int getNsteps() const {return _nsteps;}
        
        structure optimize (std::ostream &min, structure &S);
        std::pair<double,std::vector<double> > getLowestEvec (std::vector< std::pair< double,std::vector<double> > > &V);
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
                _exp2(6),
                pairPotential()
        {}

        LJ (
                double epsilon, 
                double rm, double exp1, 
                double exp2, 
                const int algo_switch, 
                const double stop_crit,
                const int nsteps        ) 
                : 
                _epsilon(epsilon),
                _rm(rm),
                _exp1(exp1),
                _exp2(exp2),
                pairPotential(algo_switch,stop_crit,nsteps)
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
        ELJ() : _c(30), pairPotential()
        {}

        ELJ (
                std::vector<double> c, 
                const int algo_switch, 
                const double stop_crit,
                const int nsteps       ) 
                : 
                _c(c),
                pairPotential(algo_switch,stop_crit,nsteps)
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
                    _range(1),
                    pairPotential()
        {}

        RangeLJ (
                    double epsilon, 
                    double rm, 
                    double exp1, 
                    double exp2, 
                    double range,
                    const int algo_switch,
                    const double stop_crit, 
                    const int nsteps       ) 
                    :   
                    _epsilon(epsilon),
                    _rm(rm),
                    _exp1(exp1),
                    _exp2(exp2),
                    _range(range),
                    pairPotential(algo_switch,stop_crit,nsteps)
        {}

        static RangeLJ *readPotential ();
};
#endif
