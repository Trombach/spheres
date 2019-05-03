#ifndef BASINHOPPING
#define BASINHOPPING

#include "structure.h"
#include "potential.h"
#include "acceptanceTest.h"

class BasinHopping
{
    public:
        BasinHopping(   structure initialCoordinates,
                        int nsteps = 1000) :  
            _initialCoordinates(initialCoordinates),
            _previousStep(initialCoordinates),
            _currentStep(initialCoordinates),
            _size(initialCoordinates.nAtoms()),
            _iteration(0),
            _nsteps(nsteps)
        {}

        int run (std::unique_ptr<AcceptanceTest>& accept);

    private:
        bool checkConf();
        bool acceptStep (double oldE, double newE, std::unique_ptr<AcceptanceTest>& accept);
        void propagate();

        const structure _initialCoordinates;
        structure _currentStep;
        structure _previousStep;
        const int _size;
        int _iteration;
        const int _nsteps;
};

#endif
