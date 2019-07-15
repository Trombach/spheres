#include "transversePotential.h"
#include "structure.h"

class EigenvectorFollowing
{
    private:
        TransversePotential _tPot;
        structure _previousStructure;
        structure _currentStructure; 
        double _max_stepsize;
        column_vector _trueGradient;
        bool _gradientIsCurrent; 

        unsigned int _currentIter;

    public:
        EigenvectorFollowing (TransversePotential &tPot, structure &S) :
            _tPot(tPot),
            _previousStructure(S),
            _currentStructure(S),
            _max_stepsize(1),
            _trueGradient(0),
            _gradientIsCurrent(false),
            _currentIter(0)
        {}

        void updateStructure (structure &S);
        column_vector getTrueGradient();
        structure getCurrentStructure() { return _currentStructure; }
        structure stepUphill (structure &S);


        void run(unsigned int n);
};
