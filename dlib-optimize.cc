#include <dlib/optimization.h>
#include <algorithm>
#include <functional>
#include "structure.h"
#include "stop_strategy.h"
#include "parameter.h"


using namespace std;
using namespace placeholders;

typedef dlib::matrix<double,0,1> column_vector;

////////////////////////////////////////////////////////////
///////////////////////////ENERGY///////////////////////////
////////////////////////////////////////////////////////////

//
//function to calculate LJ Energy, needs rij,rm and epsilon
//
double LJEnergy (const double distance, const double epsilon, const double rm, const double exp1, const double exp2) {
    return epsilon * ( (pow (rm / distance, exp1)) - 2 * (pow (rm / distance, exp2)) );
}

double structure::sumOverAllInteractions (const vector<double> &p) {
    double totalEnergy = 0;
    //iterate over double index ij, where N>j>i and N>i
    for (int i = 0; i < this->nAtoms(); ++i) { 
        for (int j = i + 1; j < this->nAtoms(); ++j) {
            totalEnergy += LJEnergy (coord3d::dist ((*this)[i],(*this)[j]), p[0], p[1], p[2], p[3]);
        }
    }
    return totalEnergy;
}




//
//function for dlib, returns f(x, params) value
//
double LJEnergy_dlib (const column_vector &v, void *params) {
    structure kissingSphere;
    for (long i = 0; i < v.size() / 3; ++i) {
        coord3d sphere(v(3 * i), v(3 * i + 1), v(3 * i + 2));
        kissingSphere.push_back(sphere);
    }
    double totalEnergy = 0;
    vector<double> *p = static_cast<vector<double>* >(params);
    for (int i = 0; i < kissingSphere.nAtoms(); ++i) { 
        for (int j = i + 1; j < kissingSphere.nAtoms(); ++j) {
            totalEnergy += LJEnergy (coord3d::dist (kissingSphere[i],kissingSphere[j]), (*p)[0], (*p)[1], (*p)[2], (*p)[3]);
        }
    }
    return totalEnergy;
}



////////////////////////////////////////////////////////////
///////////////////////////GRADIENT/////////////////////////
////////////////////////////////////////////////////////////


//
//function to calculate LJ Gradient, needs rij,rm and epsilon; returns gradient for 1 sphere pair
//
coord3d LJGradient(const coord3d ri, const coord3d rj, const double epsilon, const double rm, const double exp1, const double exp2) {
    double distance = coord3d::dist (ri , rj);
    coord3d distanceVector = rj - ri;
    double LJGradientValue = ( epsilon / rm ) * ( exp1 * (pow (rm / distance, exp1 + 1)) - 2 * exp2 * (pow (rm / distance, exp2 + 1)) );
    return distanceVector / distanceVector.norm() * LJGradientValue;
}

vector<coord3d> structure::sumOverAllGradients (const vector<double> &p) {
    vector<coord3d> gradients(this->nAtoms(), coord3d());
    coord3d force;
    for (int i = 0; i < this->nAtoms(); ++i) { 
        for (int j = i + 1; j < this->nAtoms(); ++j) {
            coord3d twoBodyGradient = LJGradient ((*this)[i], (*this)[j], p[0], p[1], p[2], p[3]);
            gradients[i] += twoBodyGradient;
            gradients[j] -= twoBodyGradient;
        }
    }
    return gradients;
}




//
//function for dlib, returns gradient of f
//
const column_vector LJGradient_dlib (const column_vector &v, void *params) {
    structure kissingSphere;
    for (long i = 0; i < v.size() / 3; ++i) {
        coord3d sphere(v(3 * i), v(3 * i + 1), v(3 * i + 2));
        kissingSphere.push_back(sphere);
    }
    vector<double> *p = static_cast<vector<double>*>(params);
    vector<coord3d> gradients = kissingSphere.sumOverAllGradients(*p);
    column_vector df(gradients.size() * 3);
    for (vector<coord3d>::size_type i = 0; i < gradients.size(); ++i) {
        for (int j=0; j <= 2; ++j) {
            df(3 * i + j) = gradients[i][j];
        }
    }
    return df;
}


////////////////////////////////////////////////////////////
/////////////////////////////OPTIMIZER/////////////////////////
////////////////////////////////////////////////////////////


//
//initialize dlib minimizer function
//
structure structure::optimize (ostream &min, parameter<int> switches, vector<double> parameters, parameter<double> opt) {
    const size_t nsteps = static_cast<const size_t>(opt.get("nsteps"));

    min.precision(16);

    const double stop_crit = opt.get("convergence");

    const int potential_switch = switches.get("potential");
    const int algo_switch = switches.get("algo");
    
    column_vector x((this->nAtoms()) * 3);
    void *params = static_cast<void*>(&parameters);
    double (*f)(const column_vector &v, void *params);
    const column_vector (*df)(const column_vector &v, void *params);
    auto f_params = bind (f, _1, params);
    auto df_params = bind (df, _1, params);

    switch (potential_switch) {
        case 1:
            f = &LJEnergy_dlib;
            df = &LJGradient_dlib;
            break;
        default:
            cerr << "Error, probably bad input of potential name!" << endl;
            structure newGeometry;
            return newGeometry;
    }


    for (int i = 0; i < this->nAtoms(); ++i) {
        for (int j = 0; j<=2; ++j) {
            x(3 * i + j) = (*this)[i][j];
        }
    }

    min << "    E                        g" << endl;

    switch (algo_switch) {
        case 1:
            try {
                dlib::find_min(dlib::bfgs_search_strategy(), 
                    dlib::stop_strategy(stop_crit, nsteps).be_verbose(min), 
                    f_params, df_params, x, -(this->nAtoms()) * 1000);
            }
            catch (std::exception &e) {
                cerr << "Structure " << this->getNumber() << ": " << e.what() << endl;
            }
            break;
        case 2: {
            //call optimizer
            try {
                dlib::find_min(dlib::cg_search_strategy(), 
                    dlib::stop_strategy(stop_crit, nsteps).be_verbose(min), 
                    f_params, df_params, x, -(this->nAtoms()) * 1000);
            }
            catch (std::exception &e) {
                cerr << "Structure " << this->getNumber() << ": " << e.what() << endl;
                structure newGeometry(this->getNumber(), this->getCoordinates());
                return newGeometry;
            }
            break;
        }
        default:
            cerr << "Error, bad input of algorithm name!" << endl;
            structure newGeometry(this->getNumber(), this->getCoordinates());
            return newGeometry;
    }
    
    vector<coord3d> coordinates;
    for (long i = 0; i < x.size() / 3; i++) {
        coord3d sphere (x(3 * i), x(3 * i + 1), x(3 * i + 2)); 
        coordinates.push_back(sphere);
    }
    
    structure newGeometry(this->getNumber(), coordinates);

    min << "-----------------------------------------------" << endl;
    min << "E: " << f_params(x) << endl;
    min << "g: " << scientific << dlib::length (df_params(x)) << endl;
    //min << "n: " << nsteps << endl;
    //if (nsteps == maxsteps) {
    //  min << "Warning, maxsteps reached!!!" << endl;
    //}
    
    newGeometry.setEnergy(f_params(x));

    return newGeometry;
}

