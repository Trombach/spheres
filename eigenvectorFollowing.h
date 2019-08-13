#include "transversePotential.h"
#include "structure.h"
#include "globals.h"

template <typename T>
class EigenvectorFollowing
{
    private:
        T _uniqueStructures;
        TransversePotential _tPot;
        structure _previousStructure;
        structure _currentStructure; 
        double _maxStepsize;
        column_vector _trueGradient;
        bool _gradientIsCurrent; 
        int _reduceStep;

        unsigned int _currentIter;

    public:
        EigenvectorFollowing (TransversePotential &tPot, structure &S) :
            _uniqueStructures(),
            _tPot(tPot),
            _previousStructure(S),
            _currentStructure(S),
            _maxStepsize(0.1),
            _trueGradient(0),
            _gradientIsCurrent(false),
            _currentIter(0),
            _reduceStep(0)
        {}

        void updateStructure (structure &S);
        column_vector getTrueGradient();
        structure getCurrentStructure() { return _currentStructure; }
        void revert() { updateStructure(_previousStructure); }
        structure stepUphill (structure &S);
        void updateMaxStepsize (structure &newS, double &overlap, double &h);

        void run(unsigned int n);
};
