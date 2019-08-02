#ifndef BASINHOPPING
#define BASINHOPPING

#include <iostream>
#include "structure.h"
#include "potential.h"
#include "acceptanceTest.h"
#include "storage.h"

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
        //StorageByEnergy _uniqueStructures;
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
