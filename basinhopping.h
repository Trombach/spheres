#ifndef BASINHOPPING
#define BASINHOPPING

#include <iostream>
#include "structure.h"
#include "potential.h"
#include "acceptanceTest.h"
#include "storage.h"

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
            _iteration(0),
            _accept(accept),
            _nsteps(nsteps),
            _uniqueStructures()
        {}

        int nStructures() {return _uniqueStructures.getSize();}
        void printEnergies(std::ostream& out) {_uniqueStructures.printEnergies(out);}

        int run ();

    private:
        bool checkConf();
        bool acceptStep (double oldE, double newE);
        void propagate();

        const structure _initialCoordinates;
        structure _currentStep;
        structure _previousStep;
        const int _size;
        int _iteration;
        const int _nsteps;
        std::shared_ptr<AcceptanceTest> _accept;
        Storage _uniqueStructures;
};

#endif
