#include <cmath>
#include "eigenvectorFollowing.h"
#include "iop.h"

using namespace std;

void EigenvectorFollowing::updateStructure (structure &S)
{
    _gradientIsCurrent = false;
    _currentStructure = S;
}

column_vector EigenvectorFollowing::getTrueGradient ()
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

structure EigenvectorFollowing::stepUphill (structure &S)
{
    double h;
    column_vector v = this->_tPot.getVector();
        double eval = fabs(_tPot.getEval());
    if (_currentIter == 0) h = 0.3;
    else
    {
        column_vector g = getTrueGradient();
        double overlap = dot(g,v);
        h = 2 * overlap / ( eval * ( 1 + sqrt( 1 + 4 * pow(2, overlap / eval) ) ) );
    }
    cout << "h: " << h << endl << "eval: " << eval << endl;

    vector<coord3d> newCoord(S.nAtoms());
    for (int i = 0; i < S.nAtoms(); i++)
        for (int j = 0; j < 3; j++)
            newCoord[i][j] = S[i][j] + h * v(3 * i + j);

    structure newS(0, newCoord, false);
    xyzout(newS, "uphill.xyz");

    return newS;
}

void EigenvectorFollowing::run (unsigned int n) 
{
        for (unsigned int i = 0; i < n; i++)
        {
            _currentIter = i;
            _tPot.calcTransverseDirection(_currentStructure);
            structure trialS = stepUphill(_currentStructure);
            structure newS = _tPot.optimize(trialS);

            _previousStructure = _currentStructure;
            updateStructure(newS);
        }
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

    ofstream dummy;
    structure S2 = potential->optimize(dummy,S1);

    vector< vector<double> > hessian = potential->calcHessian(S2);
    vector<double> eval = diag(hessian);
    for (auto& i : eval) cout << "e " << i << endl;
    xyzout(S2, "S2.xyz");
    



    //std::cout << T.getEnergy(S) << std::endl;
    //std::cout << "g2: " << T.getTrueGradient(S2) << std::endl;

    EigenvectorFollowing E(T, S2);
    E.run(1);

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

    return 0;
}
