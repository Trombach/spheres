#ifndef SPHERES
#define SPHERES

#include <iostream>
#include <vector>
#include <math.h>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multiset.h"
#include "gsl/gsl_multimin.h"


struct coord3d {
    double x[3];

    coord3d(const double y[3]) { x[0] = y[0]; x[1] = y[1]; x[2] = y[2]; }
    explicit coord3d(const double x_=0, const double y=0, const double z=0) { x[0] = x_; x[1] = y; x[2] = z; }
    coord3d operator/(const double s)   const { return coord3d(*this) /= s; }
    coord3d operator*(const double s)   const { return coord3d(*this) *= s; }
    coord3d operator*(const coord3d& y) const { return coord3d(*this) *= y; }
    coord3d operator+(const coord3d& y) const { return coord3d(*this) += y; }
    coord3d operator-(const coord3d& y) const { return coord3d(*this) -= y; }
    coord3d& operator+=(const coord3d& y){ x[0] += y[0]; x[1] += y[1]; x[2] += y[2]; return *this; }
    coord3d& operator-=(const coord3d& y){ x[0] -= y[0]; x[1] -= y[1]; x[2] -= y[2]; return *this; }
    coord3d& operator*=(const coord3d& y){ x[0] *= y[0]; x[1] *= y[1]; x[2] *= y[2]; return *this; }
    coord3d& operator*=(const double& y){ x[0] *= y; x[1] *= y; x[2] *= y; return *this; }
    coord3d& operator/=(const double& y){ x[0] /= y; x[1] /= y; x[2] /= y; return *this; }
    coord3d operator-() const {coord3d y(-x[0],-x[1],-x[2]); return y;}

	double& operator[](unsigned int i){ return x[i]; }
    double  operator[](unsigned int i) const { return x[i]; }

    double dot(const coord3d& y) const { return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; }
    double norm() const { return sqrt(dot(*this)); }
    static double dist(const coord3d& x, const coord3d& y){ return (x-y).norm(); }

    friend std::ostream& operator<<(std::ostream &s, const coord3d& x){ s << std::fixed << "{" << x[0] << "," << x[1] << "," << x[2]<< "}"; return s; }
    friend std::istream& operator>>(std::istream &s, coord3d& x){ for(int i=0;i<3;i++){ s >> x[i]; } return s; }

};


class structure:public std::vector<coord3d> {    
    
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
	//function for gsl multimin, returns f(x, params) value
	//
	double LJEnergy_gsl (const gsl_vector *v, void *params) {
		structure kissingSphere;
		for (size_t i = 0; i < v->size / 3; ++i) {
            coord3d sphere(gsl_vector_get (v, 3i), gsl_vector_get (v, 3i + 1), gsl_vector_get (v, 3i + 2));
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
            coord3d sphere(gsl_vector_get (v, 3i), gsl_vector_get (v, 3i + 1), gsl_vector_get (v, 3i + 2));
			kissingSphere.push_back(sphere);
		}
		double *p = static_cast<double*>(params);
		vector<coord3d> gradients = kissingSphere.sumOverAllGradients(p);
        for (vector<coord3d>::size_type i = 0; i < gradients.size(); ++i) {
			for (int j=0; j <= 2; ++j) {
				gsl_vector_set(df, i+j, gradients[i][j]);
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

public:   

    //
    //function to sum over all sphere interactions, change later to work with different potentials
    //
   	double sumOverAllInteractions () {
        double totalEnergy = 0;
        //iterate over double index ij, where N>j>i and N>i
    	for (structure::const_iterator iter = this->begin(); iter != this->end(); ++iter) { 
    		for (structure::const_iterator jter = iter + 1; jter != this->end(); ++jter) {
    			//cout << "iter is " << *iter << endl;
    			//cout << "jter is " << *jter << endl;
                totalEnergy += LJEnergy (coord3d::dist (*iter,*jter), 1, 0.5);
            }
    	}
    	return totalEnergy;
    }

	//
	//function to sum over all gradients to get gradient for each sphere
	//
	vector<coord3d> sumOverAllGradients (const double *p) {
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
	//initialize gsl minimizer function
	//
	void optimize () {
		gsl_multimin_function_fdf min_function;

		double p[2] = { 1.0, 0.5 };

		min_function.n = (this->size()) * 3;
        min_function.f = &structure::LJEnergy_gsl;
		min_function.df = &LJGradient_gsl;
		min_function.fdf = &LJEnergyAndGradient_gsl;
		min_function.params = static_cast<void*>(p);
	}
};

#endif
