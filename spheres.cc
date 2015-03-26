#include "gsl/gsl_vector.h"
#include "gsl/gsl_multimin.h"
#include "spheres.h"


using namespace std;

//
//function to calculate LJ Energy, needs rij,rm and epsilon
//
double LJEnergy (const double distance, const double epsilon, const double rm) {
    return epsilon * ( (pow (rm / distance, 12)) - 2 * (pow (rm / distance, 6)) );
}


//
//function to calculate LJ Gradient, needs rij,rm and epsilon; returns gradient for 1 sphere pair
//
coord3d LJGradient(const coord3d ri, const coord3d rj, const double epsilon, const double rm) {
	double distance = coord3d::dist (ri , rj);
	coord3d distanceVector = rj - ri;
	double LJGradientValue = ( 12 * epsilon / rm ) * ( (pow (rm / distance, 13)) - (pow (rm / distance, 7)) );
	return distanceVector / distanceVector.norm() * LJGradientValue;
}

//
//function to calculate LJ Hessian
//
double LJHessian (const double distance, const double epsilon, const double rm) {
    return - 12 * epsilon / pow (rm, 2) * (13 * pow (rm / distance, 14) - 7 * pow (rm / distance, 8) );
}

double structure::sumOverAllInteractions () {
    double totalEnergy = 0;
    //iterate over double index ij, where N>j>i and N>i
	for (structure::const_iterator iter = this->begin(); iter != this->end(); ++iter) { 
		for (structure::const_iterator jter = iter + 1; jter != this->end(); ++jter) {
			//cout << "iter is " << *iter << endl;
			//cout << "jter is " << *jter << endl;
            totalEnergy += LJEnergy (coord3d::dist (*iter,*jter), 1.0, 1.0);
        }
	}
	return totalEnergy;
}


vector<coord3d> structure::sumOverAllGradients (const double *p) {
	vector<coord3d> gradients(this->size(), coord3d());
	coord3d force;
	for (structure::size_type i = 0; i < this->size(); ++i) { 
		for (structure::size_type j = i + 1; j < this->size(); ++j) {
		    coord3d twoBodyGradient = LJGradient ((*this)[i], (*this)[j], p[0], p[1]);
            gradients[i] += twoBodyGradient;
			gradients[j] -= twoBodyGradient;
        }
	}
    return gradients;
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
	double *p = static_cast<double*>(params);
	for (structure::const_iterator iter = kissingSphere.begin(); iter != kissingSphere.end(); ++iter) { 
		for (structure::const_iterator jter = iter + 1; jter != kissingSphere.end(); ++jter) {
			//cout << "iter is " << *iter << endl;
			//cout << "jter is " << *jter << endl;
            totalEnergy += LJEnergy (coord3d::dist (*iter,*jter), p[0], p[1]);
        }
	}
	return totalEnergy;
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
	double *p = static_cast<double*>(params);
	vector<coord3d> gradients = kissingSphere.sumOverAllGradients(p);
    for (vector<coord3d>::size_type i = 0; i < gradients.size(); ++i) {
		for (int j=0; j <= 2; ++j) {
			gsl_vector_set(df, 3 * i + j, gradients[i][j]);
		}
	}
}


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
void structure::optimize (vector<double> parameters, const vector<double> opt) {
	int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
    
	gsl_vector *x;
	gsl_multimin_function_fdf min_function;

    double p[2] = { 1.0, 1.0 };

	min_function.n = (this->size()) * 3;
    min_function.f = &LJEnergy_gsl;
	min_function.df = &LJGradient_gsl;
	min_function.fdf = &LJEnergyAndGradient_gsl;
	min_function.params = static_cast<void*>(p);

	x = gsl_vector_alloc ( (this->size()) *3 );
    for (structure::size_type i = 0; i < this->size(); ++i) {
	    for (int j = 0; j<=2; ++j) {
	        gsl_vector_set(x, 3 * i + j, (*this)[i][j]);
		}
	}

	T = gsl_multimin_fdfminimizer_vector_bfgs2;
	s = gsl_multimin_fdfminimizer_alloc (T, (this->size()) * 3);
    
	double stepSize = opt[2];
	double accuracy = opt[0];
	gsl_multimin_fdfminimizer_set (s, &min_function, x, stepSize, accuracy);

	size_t i = 0;
	do {
		cout << "-----------------------------------------------" << endl;
		i++;
		status = gsl_multimin_fdfminimizer_iterate (s);
        cout << status << endl;
		
		if (status)
			break;
        
		double absoluteTolerance = opt[1];
		status = gsl_multimin_test_gradient (s->gradient, absoluteTolerance);
        //cout << status << endl;

		if (status == GSL_SUCCESS)
			cout << "Minimum found at:\n" << endl;

		//create structure for optimized geometry
    	structure kissingSphere;
    	for (size_t i = 0; i < x->size / 3; ++i) {
            coord3d sphere(gsl_vector_get (s->x, 3 * i), gsl_vector_get (s->x, 3 * i + 1), gsl_vector_get (s->x, 3 * i + 2));
			cout << sphere <<  endl;
    		kissingSphere.push_back(sphere);
    	}
    cout << "LJ Energy is: " << s->f << endl;
	}
	while (status == GSL_CONTINUE && i < 100);

    gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);

}
