#include <dlib/optimization.h>
#include <cassert>
#include "stop_strategy.h"
#include "potential.h"
#include "transversePotential.h"
#include "lina.h"
#include "iop.h"

using namespace std;

void TransversePotential::setVector (double &eval, column_vector &v) 
{ 
    _eval = eval;
    _vector = v / sqrt(dot(v, v)); 
}

void TransversePotential::setVector (double &eval, vector<double> &v)
{
    _eval = eval;
    column_vector vector(v.size());
    for (int i = 0; i < v.size() ; i++) 
    {
        vector(i) = v[i];
    }

    _vector = vector / sqrt(dot(vector, vector));
}

column_vector TransversePotential::getTransverseGradient (const column_vector &v)
{
    assert(_vector.size() == v.size());
    const column_vector g = getTrueGradient(v);
    assert(_vector.size() == g.size());

    return g - dot(g, _vector) * _vector; 

}

column_vector TransversePotential::getTransverseGradient (structure &S)
{
    assert(_vector.size()); 
    column_vector g = getTrueGradient(S);
    assert(_vector.size() == g.size());

    return g - dot(g, _vector) * _vector;
}

void TransversePotential::calcTransverseDirection( structure &S )
{
    vector< vector<double> > hessian = _potential->calcHessian(S);
    
    vector<pair<double, vector<double> > > eval_evec = diagv(hessian);

    pair<double,vector<double> > lowest_eval_evec = _potential->getLowestEvec(eval_evec);
    setVector(lowest_eval_evec.first, lowest_eval_evec.second);
}


structure TransversePotential::optimize (structure &S)
{
    //calcTransverseDirection(S);

    column_vector x((S.nAtoms()) * 3);
    for (int i = 0; i < S.nAtoms(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            x(3 * i + j) = S[i][j];
        }
    }
    

    auto f = [this] (column_vector v) -> const double {return this->getEnergy(v);};
    auto df = [this] (column_vector v) -> column_vector {return this->getTransverseGradient(v);};

    try
    {
        switch (_potential->getAlgoSwitch())
        {
            case 1:
            {
                dlib::find_min( dlib::bfgs_search_strategy(),
                                dlib::stop_strategy(_potential->getStopCrit(), _potential->getNsteps()).be_verbose(),
                                f, df, x, -(S.nAtoms()) * 1000);
                break;
            }
            case 2:
            {
                dlib::find_min( dlib::cg_search_strategy(),
                                dlib::stop_strategy(_potential->getStopCrit(), _potential->getNsteps()),
                                f, df, x, -(S.nAtoms()) * 1000);
                break;
            }
            default:
            {
                cerr << "Bad input of algorithm name" << endl;
                structure newS(S.getNumber(), S.getCoordinates());
                return newS;
            }
        }
    }
    catch (std::exception &e)
    {
        structure newS(S.getNumber(), S.getCoordinates());
        return newS;
    }

    vector<coord3d> newCoordinates;
    for (long i = 0; i < x.size() / 3; i++)
    {
        coord3d sphere (x(3 * i), x(3 * i + 1), x(3 * i + 2));
        newCoordinates.push_back(sphere);
    }

    structure newS(S.getNumber(), newCoordinates);
    double finalEnergy = this->getEnergy(x);

    newS.setEnergy(finalEnergy);
    newS.setConverged();

    cout << "gT: " << dlib::length(this->getTransverseGradient(x)) << endl;
    cout << "g: " << dlib::length(this->getTrueGradient(x)) << endl;
    vector< vector<double> > hessian = _potential->calcHessian(newS);
    vector<double> eval = diag(hessian);
    for (auto& i : eval)
        cout << "H: " << i << endl;

    return newS;
}

