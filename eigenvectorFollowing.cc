#include <cmath>
#include "eigenvectorFollowing.h"
#include "iop.h"
#include "globals.h"
#include "storage.h"

using namespace std;

template <typename T>
void EigenvectorFollowing<T>::updateStructure (structure &S)
{
    _gradientIsCurrent = false;

    vector< vector<double> > hessian = _tPot.getTrueHessian(S);
    vector<double> eval = diag(hessian);
    S.setHessian(eval);

    _currentStructure = S;
}

template <typename T>
column_vector EigenvectorFollowing<T>::getTrueGradient ()
{
    if (_gradientIsCurrent) 
        return _trueGradient;
    else
    {
        _gradientIsCurrent = true;
        _trueGradient = _tPot.getTrueGradient(_currentStructure);
        return _trueGradient;
    }
}

template <typename T>
void EigenvectorFollowing<T>::updateMaxStepsize (structure &newS, double &overlap, double &h) 
{
    column_vector newGrad = _tPot.getTrueGradient(newS);
    double eval = _tPot.getEval();
    double newOverlap = dot(newGrad, _tPot.getVector());
    double a = 1 - (newOverlap - overlap) / (h * eval);
    double b = 1 - (- newOverlap - overlap) / (h * eval);
    double val = min(fabs(a),fabs(b));

    if (val > 2) _maxStepsize = max(_maxStepsize / 1.1, 0.01);
    else _maxStepsize = min(_maxStepsize * 1.1, 0.5);

    cout << "max stepsize: " << _maxStepsize << endl;
}

template <typename T>
structure EigenvectorFollowing<T>::stepUphill (structure &S)
{
    double h;
    double eval = _tPot.getEval();
    column_vector v = this->_tPot.getVector();
    column_vector g = getTrueGradient();
    double overlap = dot(g,v);
    double maxStepsize = _maxStepsize;

    if (_reduceStep > 0) maxStepsize *= pow(0.1,_reduceStep);

    if (_currentIter == 0 || eval > 0.001) h = maxStepsize;
    else 
    {
        double absval = fabs(eval);
        h = 2 * overlap / ( absval * ( 1 + sqrt( 1 + 4 * pow(2, overlap / absval) ) ) );
    }
    cout << "h: " << h << endl;
    cout << "eval: " << _tPot.getEval() << endl;
    //cout << "v: " << _tPot.getVector() << endl;
    cout << length(_tPot.getVector()) << endl;
    vector< vector<double> > hessian = _tPot.getTrueHessian( S);
    vector<double> evalhess = diag(hessian);
    for (auto& i : evalhess) cout << "H: " << i << endl; 
    //cout << "length: " << v << endl;
    //cout << "length: " << g << endl;
    
    if (fabs(h) > maxStepsize) h *= maxStepsize / fabs(h);

    vector<coord3d> newCoord(S.nAtoms());
    for (int i = 0; i < S.nAtoms(); i++)
        for (int j = 0; j < 3; j++)
            newCoord[i][j] = S[i][j] + h * v(3 * i + j);


    structure newS(0, newCoord, false);
    if (eval < 0) updateMaxStepsize(newS, overlap, h);
    xyzout(newS, "uphill.xyz");

    return newS;
}

template <typename T>
void EigenvectorFollowing<T>::run (unsigned int n) 
{
        for (unsigned int i = 0; i < n; i++)
        {
            cout << "------------------" << endl << "Step " << i << endl;
            _currentIter = i;
            _tPot.calcTransverseDirection(_currentStructure);

            //checkup
            if (_tPot.getEval() > 0 && _currentIter > 0)
            {
                cerr << "reducing step!" << endl;
                cerr << _tPot.getEval() << " " << _currentIter << endl;
                revert();
                _reduceStep++;
                _tPot.calcTransverseDirection(_currentStructure);
            }
            else _reduceStep = 0;

            cout << "Reduce step: " << _reduceStep << endl;

            structure trialS = stepUphill(_currentStructure);
            double currentE = _tPot.getEnergy(_currentStructure);
            double uphillE = _tPot.getEnergy(trialS);
            if ( uphillE < currentE ) 
            {
                cerr << "warning: uphill step decreased energy, " << currentE << " -> " << uphillE << endl;
            }
            cout << "current E: " << currentE << endl;
            cout << "optimizing..." << endl;
            structure newS = _tPot.optimize(trialS);
            cout << "new E: " <<_tPot.getEnergy(newS) << endl;
            

            _previousStructure = _currentStructure;
            updateStructure(newS);

            _uniqueStructures.addCluster(newS);

            if (_currentStructure.countNegativeEval() < 1) continue;
            if (fabs(length(_tPot.getTrueGradient(newS))) < 1e-5) break;
        }

        _uniqueStructures.printStructures();
        cout << _uniqueStructures.getSize() << endl;
        _uniqueStructures.printKeys(cout);
        cout <<"done"<< endl;
}

//------------------------------------------------------------------------------//
//                           main for testing                                   //
//------------------------------------------------------------------------------//


int main ()
{
    
    LJ *potential = new LJ();
    TransversePotential T(*potential);

    vector<coord3d> coords;
    /*
    coords.push_back(coord3d(0.25877219650832,-0.51611072802221,0));
    coords.push_back(coord3d(0.31757890372974,0.48215865980443,0));
    coords.push_back(coord3d(-0.57635110023806,0.033952068217783,0));
    */
    coords.push_back(coord3d(1,0,0));
    coords.push_back(coord3d(-1,0,0));
    coords.push_back(coord3d(0,1,0));
    coords.push_back(coord3d(0,-1,0));
    coords.push_back(coord3d(0,0,1));
    coords.push_back(coord3d(0,0,-1));
    structure S1(1,coords);

    vector<coord3d> coord2;
    coord2.push_back(coord3d(0,0,0));
    coord2.push_back(coord3d(0,0,1));
    structure lin(2, coord2);

    coord2.push_back(coord3d(0.5,0.5,0));
    structure tri(3, coord2);

    coord2.push_back(coord3d(0.5,0.25,0.5));
    //coord2.push_back(coord3d(0.5,0.25,-0.5));
    structure bipyr(4, coord2);

    structure rand(6);

    ofstream dummy;
    structure S2 = potential->optimize(dummy, rand);

    vector< vector<double> > hessian = potential->calcHessian(S2);
    vector<double> eval = diag(hessian);
    for (auto& i : eval) cout << "e " << i << endl;
    xyzout(S2, "S2.xyz");
    



    //std::cout << T.getEnergy(S) << std::endl;
    //std::cout << "g2: " << T.getTrueGradient(S2) << std::endl;

    EigenvectorFollowing<StorageByInterPartDist> E(T, S2);
    E.run(4);

    structure newS = E.getCurrentStructure();
    xyzout(newS, "newS.xyz");
    
    /*
    vector<coord3d> coordinates = newS.getCoordinates();
    for (auto& i : coordinates)
        cout << i << endl;
    */

    /*
    column_vector gradient = potential->calcGradient(newS);
    column_vector gradientT = T.getTransverseGradient(newS);
    vector< vector<double> > hessian = potential->calcHessian(newS);
    vector<double> eval = diag(hessian);
    */

    //xyzout(newS, "newS.xyz");

    /*
    for (auto& i : gradientT)
        cout << "gT: " << i << endl;
    for (auto& i : gradient)
        cout << "g: " << i << endl;
    for (auto& i : eval)
        cout << "H2: " << i << endl;
    */
    delete potential;

    return 0;
}
