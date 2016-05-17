#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <algorithm>
#include "structure.h"


using namespace std;

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
//function for gsl multimin, returns f(x, params) value
//
double LJEnergy_gsl (const gsl_vector *v, void *params) {
	structure kissingSphere;
	for (size_t i = 0; i < v->size / 3; ++i) {
        coord3d sphere(gsl_vector_get (v, 3 * i), gsl_vector_get (v, 3 * i + 1), gsl_vector_get (v, 3 * i + 2));
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
//function for gsl multimin, returns gradient of f
//
void LJGradient_gsl (const gsl_vector *v, void *params, gsl_vector *df) {
	structure kissingSphere;
	for (size_t i = 0; i < v->size / 3; ++i) {
        coord3d sphere(gsl_vector_get (v, 3 * i), gsl_vector_get (v, 3 * i + 1), gsl_vector_get (v, 3 * i + 2));
		kissingSphere.push_back(sphere);
	}
	vector<double> *p = static_cast<vector<double>*>(params);
	vector<coord3d> gradients = kissingSphere.sumOverAllGradients(*p);
    for (vector<coord3d>::size_type i = 0; i < gradients.size(); ++i) {
		for (int j=0; j <= 2; ++j) {
			gsl_vector_set(df, 3 * i + j, gradients[i][j]);
		}
	}
}


////////////////////////////////////////////////////////////
/////////////////////////////OPTIMIZER/////////////////////////
////////////////////////////////////////////////////////////



//
//function for gsl multimin, compute f and df together
//
void LJEnergyAndGradient_gsl (const gsl_vector *x, void *params, double *f, gsl_vector *df) {
	*f = LJEnergy_gsl(x, params);
	LJGradient_gsl(x, params, df);
}


//
//initialize gsl minimizer function
//
structure structure::optimize (ofstream &min, const int &algo_switch, const int &potential_switch, vector<double> parameters, const vector<double> opt, vector<double> &allEnergies) {
	int status;
	structure newGeometry;
	size_t nsteps = static_cast<size_t>(opt[3]);

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
    
	gsl_vector *x;
	gsl_multimin_function_fdf min_function;
    

    switch (potential_switch) {
		case 1:
            min_function.f = &LJEnergy_gsl;
	        min_function.df = &LJGradient_gsl;
	        min_function.fdf = &LJEnergyAndGradient_gsl;
			break;
		default:
			cerr << "Error, probably bad input of potential name!" << endl;
			return newGeometry;
    }

    min_function.n = (this->nAtoms()) * 3;
	min_function.params = static_cast<void*>(&parameters);

	x = gsl_vector_alloc ( (this->nAtoms()) *3 );
    for (int i = 0; i < this->nAtoms(); ++i) {
	    for (int j = 0; j<=2; ++j) {
	        gsl_vector_set(x, 3 * i + j, (*this)[i][j]);
		}
	}


	switch (algo_switch) {
		case 1:
	        T = gsl_multimin_fdfminimizer_vector_bfgs2;
			break;
		default:
			cerr << "Error, bad input of algorithm name!" << endl;
			return newGeometry;
	}

	s = gsl_multimin_fdfminimizer_alloc (T, (this->nAtoms()) * 3);
    
	double stepSize = opt[2];
	double accuracy = opt[0];
	gsl_multimin_fdfminimizer_set (s, &min_function, x, stepSize, accuracy);

	size_t i = 0;
	do {
		i++;
		status = gsl_multimin_fdfminimizer_iterate (s);
		
		if (status) {
			cerr << "GSL error:\n\t" << status << ": " << gsl_strerror(status) << endl;
			break;
		}
        
		double absoluteTolerance = opt[1];
		status = gsl_multimin_test_gradient (s->gradient, absoluteTolerance);

		if (status == GSL_SUCCESS) {
			min << "Minimum found at:\n" << endl;

		    //create structure for optimized geometry
    	    for (size_t j = 0; j < x->size / 3; ++j) {
                coord3d sphere(gsl_vector_get (s->x, 3 * j), gsl_vector_get (s->x, 3 * j + 1), gsl_vector_get (s->x, 3 * j + 2));
			    min << sphere <<  endl;
    		    newGeometry.push_back(sphere);
    	    }
        allEnergies.push_back(s->f);
		}
	}
	while (status == GSL_CONTINUE && i <= nsteps);

    
    min << "-----------------------------------------------" << endl;
    min << "LJ Energy is: " << s->f << endl;
	newGeometry.setEnergy(s->f);

    if (status == GSL_SUCCESS) {
		min << "Optimization successful!" << endl;
	}
	else {
		min << "Error: " << status << ". Optimization failed." << endl;
	}
    gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);

    return newGeometry;
}

