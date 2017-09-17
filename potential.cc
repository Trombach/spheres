#include <dlib/optimization.h>
#include <vector>
#include <iostream>
#include <libconfig.h++>
#include "geometry.h"
#include "structure.h"
#include "parameter.h"
#include "stop_strategy.h"
#include "potential.h"
#include "iop.h"

using namespace std;


/*--------------------------------------------------------------------------------------*/
//                          pair potential base class
/*--------------------------------------------------------------------------------------*/

pairPotential::pairPotential (void) {}

/*----------------------------------Energy----------------------------------------------*/

double pairPotential::calcEnergy (const column_vector &v)
{
    structure S;
    for (long i = 0; i < v.size() / 3; ++i) 
    {
        coord3d sphere(v(3 * i), v(3 * i + 1), v(3 * i + 2));
        S.push_back(sphere);
    }
    double f(0);
    for (int i = 0; i < S.nAtoms(); ++i) {
        for (int j = i + 1; j < S.nAtoms(); ++j) {
            f += this->E (coord3d::dist (S[i],S[j]));
        }
    }
    return f;
}

/*----------------------------------Energy from structure-------------------------------*/

double pairPotential::calcEnergy (structure &S)
{
    double f(0);
    for (int i = 0; i < S.nAtoms(); ++i) {
        for (int j = i + 1; j < S.nAtoms(); ++j) {
            f += this->E (coord3d::dist (S[i],S[j]));
        }
    }
    return f;
}


/*----------------------------------Gradient--------------------------------------------*/

const column_vector pairPotential::calcGradient (const column_vector &v)
{
    structure S;
    for (long i = 0; i < v.size() / 3; ++i) 
    {
        coord3d sphere(v(3 * i), v(3 * i + 1), v(3 * i + 2));
        S.push_back(sphere);
    }
    vector<coord3d> gradients (S.nAtoms(), coord3d());
    for (int i = 0; i < S.nAtoms(); i++)
    {
        for (int j = i + 1; j < S.nAtoms(); j++)
        {
            coord3d distanceVector = S[i]-S[j];
            double gradValue = this->dE_dr (coord3d::dist(S[i],S[j]));
            coord3d twoBodyGradient = distanceVector / distanceVector.norm() * gradValue;
            gradients[i] += twoBodyGradient;
            gradients[j] -= twoBodyGradient;
        }
    }
    column_vector df (gradients.size() * 3);
    for (vector<coord3d>::size_type i = 0; i< gradients.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            df(3 * i + j) = gradients[i][j];
        }
    }
    return df;
}

/*----------------------------------Gradient from structure-----------------------------*/

const column_vector pairPotential::calcGradient (structure &S)
{
    vector<coord3d> gradients (S.nAtoms(), coord3d());
    for (int i = 0; i < S.nAtoms(); i++)
    {
        for (int j = i + 1; j < S.nAtoms(); j++)
        {
            coord3d distanceVector = S[i]-S[j];
            double gradValue = this->dE_dr (coord3d::dist(S[i],S[j]));
            coord3d twoBodyGradient = distanceVector / distanceVector.norm() * gradValue;
            gradients[i] += twoBodyGradient;
            gradients[j] -= twoBodyGradient;
        }
    }
    column_vector df (gradients.size() * 3);
    for (vector<coord3d>::size_type i = 0; i< gradients.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            df(3 * i + j) = gradients[i][j];
        }
    }
    return df;
}

/*----------------------------------Hessian---------------------------------------------*/

vector< vector<double> > pairPotential::calcHessian (structure &S)
{
    vector< vector<double> > hessianMatrix (S.nAtoms() * 3, vector<double> (S.nAtoms() * 3, 0));

    for (int i = 0; i < S.nAtoms(); i++)
    {
        for (int j = i + 1; j < S.nAtoms(); j++)
        {
            const coord3d vecr = S[i] - S[j];
            const double r = coord3d::dist(S[i], S[j]);

            //calculate first and second derivative values
            double dE_dr = this->dE_dr(r);
            double d2E_dr2 = this->d2E_dr2(r);

            //calculate derivatives of r
            coord3d dvecr_dr = coord3d::dnorm(vecr);
            vector<double> d2rvecr_dr2(9, double());
            coord3d::ddnorm(vecr, d2rvecr_dr2);

            //calculation of hessian elements
            //loop over all 6 coordinates of 1 atom pair
            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    //calculate the value first, which will always only differ by sign
                    const double hessianValue = dE_dr * d2rvecr_dr2[3 * k + l] 
                                                + d2E_dr2 * dvecr_dr[k] * dvecr_dr[l];

                    /*write hessian
                     this is basically a 2 atom hessian, where the diagonal quadrants are the same and the
                     remaining quadrants are of the opposite sign
                     each quadrant can have contributions from different atom pairs
                     eg atom pair 1/2 and 1/3 both have non zero second derivatives with respect to the 
                     coordinates of atom 1*/
                    hessianMatrix[3 * i + k][3 * i + l] += hessianValue;
                    hessianMatrix[3 * i + k][3 * j + l] -= hessianValue;
                    hessianMatrix[3 * j + k][3 * i + l] -= hessianValue;
                    hessianMatrix[3 * j + k][3 * j + l] += hessianValue;
                }
            }
        }
    }
    return hessianMatrix;
}


/*----------------------------------Optimization----------------------------------------*/

structure pairPotential::optimize (ostream &min, structure &S, parameter<int> &switches, parameter<double> &opt)
{
    min.precision(16);

    const size_t nsteps = static_cast<const size_t>(opt.get("nsteps"));
    const double stop_crit = opt.get("convergence");

    const int algo_switch = switches.get("algo");

    column_vector x((S.nAtoms()) * 3);
    for (int i = 0; i < S.nAtoms(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            x(3 * i + j) = S[i][j];
        }
    }

    
    auto f = [this] (column_vector v) -> const double { return this->calcEnergy(v); };
    auto df = [this] (column_vector v) -> column_vector { return this->calcGradient(v); };


    switch (algo_switch)
    {
        case 1:
            {
                try
                {
                    dlib::find_min( dlib::bfgs_search_strategy(),
                                    dlib::stop_strategy(stop_crit, nsteps).be_verbose(min),
                                    f, df, x, -(S.nAtoms()) * 1000);
                }
                catch (std::exception &e)
                {
                    cerr << "Structure " << S.getNumber() << ": " << e.what() << endl;
                    structure newS(S.getNumber(), S.getCoordinates());
                    return newS;
                }
                break;
            }
        case 2:
            {
                try
                {
                    dlib::find_min( dlib::cg_search_strategy(),
                                    dlib::stop_strategy(stop_crit, nsteps).be_verbose(min),
                                    f, df, x, -(S.nAtoms()) * 1000);
                }
                catch (std::exception &e)
                {
                    cerr << "Structure " << S.getNumber() << ": " << e.what() << endl;
                    structure newS(S.getNumber(), S.getCoordinates());
                    return newS;
                }
                break;
            }
        default:
            {
                cerr << "Bad input of algorithm name" << endl;
                structure newS(S.getNumber(), S.getCoordinates());
                return newS;
            }
    }

    vector<coord3d> newCoordinates;
    for (long i = 0; i < x.size() / 3; i++)
    {
        coord3d sphere (x(3 * i), x(3 * i + 1), x(3 * i + 2));
        newCoordinates.push_back(sphere);
    }

    structure newS(S.getNumber(), newCoordinates);
    double finalEnergy = this->calcEnergy(x);
    column_vector finalGradient = this->calcGradient(x);

    min << "-----------------------------------------------" << endl;
    min << "E: " << finalEnergy << endl;
    min << "g: " << scientific << dlib::length(finalGradient) << endl;

    newS.setEnergy(finalEnergy);


    return newS;
}


/*--------------------------------------------------------------------------------------*/
//                          LJ potential derived class
/*--------------------------------------------------------------------------------------*/


double LJ::E (double distance)
{
    return (_epsilon / (_exp1/_exp2 - 1)) * ( (pow (_rm / distance, _exp1)) - (_exp1/_exp2) * (pow (_rm / distance, _exp2)) );
}

double LJ::dE_dr (double distance)
{
    return - (( _epsilon / (_rm * (_exp1/_exp2 - 1)) ) * ( _exp1 * (pow (_rm / distance, _exp1 + 1)) - _exp1 * (pow (_rm / distance, _exp2 + 1)) ));
}

double LJ::d2E_dr2 (double distance)
{
    return _epsilon / (pow (_rm, 2) * (_exp1/_exp2-1)) * ( (pow (_exp1, 2) + _exp1) * pow (_rm / distance, _exp1 + 2) - (_exp1 * _exp2 + _exp1) * pow (_rm / distance, _exp2 + 2) );
}

LJ *LJ::readPotential ()
{
    libconfig::Config cfg;
    try {
        cfg.readFile("settings");
    }
    catch(const libconfig::FileIOException &fioex) {
        cerr << "\tI/O error while reading file." << endl;
    }
    catch(const libconfig::ParseException &pex) {
        cerr << "\tParse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << endl;
    }
    double epsilon, rm, exp1, exp2;
    cfg.lookupValue("potential.epsilon", epsilon);
    cfg.lookupValue("potential.rm", rm);
    cfg.lookupValue("potential.exp1", exp1);
    cfg.lookupValue("potential.exp2", exp2);

    LJ *potential = new LJ(epsilon, rm, exp1, exp2); 

    return potential;
}


/*--------------------------------------------------------------------------------------*/
//                          extended LJ potential derived class
/*--------------------------------------------------------------------------------------*/

double ELJ::E (double distance)
{
    double en(0);
    for (vector<double>::size_type i = 0; i < _c.size(); i++)
    {
        en += _c[i] * pow(1 / distance, i);
    }
    return en;
}

double ELJ::dE_dr (double distance)
{
    double grad(0);
    for (vector<double>::size_type i = 0; i < _c.size(); i++)
    {
        grad -=  i * _c[i] * pow(1 / distance, i + 1);
    }
    return grad;
}

double ELJ::d2E_dr2 (double distance)
{
    double hess(0);
    for (vector<double>::size_type i = 0; i < _c.size(); i++)
    {
        hess += i * (i + 1) * _c[i] * pow(1 / distance, i + 2);
    }
    return hess;
    }

ELJ *ELJ::readPotential ()
{
    ifstream input;
    vector<double> c(30);

    if (!fexists("ext"))
    {
        cerr << "No extended LJ parameters found" << endl;
    }
    input.open("ext");

    int n;
    double c_value;
    while (input >> n >> c_value)
    {
        c[n] = c_value;
    }

    ELJ *potential = new ELJ(c);

    return potential;
}


/*--------------------------------------------------------------------------------------*/
//                          main function for testing
/*--------------------------------------------------------------------------------------*/


//int main ()
//{
//    vector<coord3d> coords;
//    coords.push_back(coord3d(0,0,0));
//    coords.push_back(coord3d(1,0,0));
//    coords.push_back(coord3d(1,0.87,0));
//    structure test(1,coords);
//
//    column_vector x(test.nAtoms() * 3);
//    for (int i = 0; i < test.nAtoms(); ++i) 
//    {
//        for (int j = 0; j<=2; ++j) 
//        {
//            x(3 * i + j) = (test)[i][j];
//        }
//    }
//    pairPotential pp;
//    LJ len;
//    cout << len.calcEnergy(x) << endl;
//    const column_vector grad = len.calcGradient(x);
//    cout << grad << endl;
//
//
//    //len.optimize<LJ>(std::cout, test);
//    //pp.optimize<pairPotential>(std::cout, test);
//    
//    len.calcHessian(test);
//
//return 0;
//}

