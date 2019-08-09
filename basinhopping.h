#ifndef BASINHOPPING
#define BASINHOPPING

#include <iostream>
#include "structure.h"
#include "potential.h"
#include "acceptanceTest.h"
#include "storage.h"

/*
This class implements monte-carlo basin hopping and stores unique local minima.

Implementations can be found in "basinhopping.cc".

Parameters
----------

const structure _initialCoordinates :
    saves coordinates of initial structure
structure _currentStep :
    saves coordinates of current step
structure _previousStep :
    saves coordinates of previous step
const int _size :
    equal to the number of atoms in the initialising structure
const int _nsteps :
    maximum number of basinhopping steps
std::shared_ptr<AcceptanceTest> _accept :
    acceptance test. See "acceptanceTest.h".
T _uniqueStructures :
    storage for unique structures. See "storage.h".
bool _accepted :
    true if Monte-Carlo step was accepted, false if otherwise.
unsigned int _iteration :
    current number of iterations.
unsigned int _naccept :
    number of accepted configurations (structures).
unsigned int _nattempts :
    number of times the stepsize was attempted to be updated.
unsigned int _nsame :
    number of times the same structure has been found.
int _interval :
    number of steps that defines how often temperature and stepsize should be adjusted.
double _stepScale :
    scaling parameter for the stepsize. dynamically adjusted.
*/

template <typename T>
class BasinHopping
{
    public:
        BasinHopping(   structure initialCoordinates,
                        std::shared_ptr<AcceptanceTest>& accept, 
                        int nsteps = 1000) :  
            _initialCoordinates(initialCoordinates),
            _previousStep(initialCoordinates),
            _currentStep(initialCoordinates),
            _size(initialCoordinates.nAtoms()),
            _accept(accept),
            _nsteps(nsteps),
            _uniqueStructures(),
            _accepted(false),
            _iteration(0),
            _naccept(0),
            _nattempts(0),
            _nsame(0),
            _interval(300),
            _stepScale(0.4)
        {}

        int nStructures() {return _uniqueStructures.getSize();}
        void printEnergies(std::ostream& out) {_uniqueStructures.printKeys(out);}

        int run ();

    private:
        bool checkConf();
        bool acceptStep (double oldE, double newE);
        void propagate();

        const structure _initialCoordinates;
        structure _currentStep;
        structure _previousStep;
        const int _size;
        const int _nsteps;
        std::shared_ptr<AcceptanceTest> _accept;
        T _uniqueStructures;

        bool _accepted;
        unsigned int _iteration;
        unsigned int _naccept;
        unsigned int _nattempts;
        unsigned int _nsame;
        int _interval;
        double _stepScale;

        void updateStep();
        void adjustStep();
        void adjustTemp();
        void resetUpdateStep();
};
#endif
