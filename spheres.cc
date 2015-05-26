#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
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


vector< vector<double> > structure::hessian (const vector<double> &p) {
	const double epsilon = p[0];
	const double rm = p[1];
	vector < vector<double>  > hessianMatrix (this->size() * 3, vector<double> (this->size() * 3, 0));
    for (structure::size_type i = 0; i < this->size(); i++) {
		for (structure::size_type j = i+1; j < this->size(); j++) {

			const coord3d vecr = (*this)[i] - (*this)[j];
			const double r = coord3d::dist((*this)[i], (*this)[j]);

			//calculate first derivative
			double dE_dr = ( 12 * epsilon / rm ) * ( (pow (rm / r, 13)) - (pow (rm / r, 7)) );

			//calculate second derivative
            double d2E_dr2 = 12 * epsilon / pow (rm, 2) * (13 * pow (rm / r, 14) - 7 * pow (rm / r, 8) );

			//calculate derivatives of r
			coord3d dvecr_dr = coord3d::dnorm(vecr);
			//cout << "dvecr_dr" << endl;
			//for (int i = 0; i < 3; i++) {cout << dvecr_dr[i] << " ";}
			//cout << endl << "------------------" <<endl;

			vector<double> d2rvecr_dr2(9, double());
			coord3d::ddnorm(vecr, d2rvecr_dr2);


			//calculate dr_dx terms
			//dnorm gives dvecr_dr which can be transformed into dr_dx by
			//   dr_dx = dr_dri * dri_dx, where dri is one element of vecr
			//the last term gives either +1 or -1
			//vector<double> dr_dx;
			//for (int a = 0; a < 3; a++) {
			//    dr_dx.push_back(dvecr_dr[a]);
			//}
			//for (int a = 0; a < 3; a++) {
			//    dr_dx.push_back(- dvecr_dr[a]);
			//}

			//calculate d2r_dx2 terms
			
			//vector<double> d2r_dx2;
			//for (int a = 0; a < 9; a++) {
			//	d2r_dx2.push_back(d2rvecr_dr2[a]);
			//}
			//for (int a = 0; a < 9; a++) {
			//	d2r_dx2.push_back(- d2rvecr_dr2[a]);
			//}
			//for (int a = 0; a < 9; a++) {
			//	d2r_dx2.push_back(- d2rvecr_dr2[a]);
			//}
			//for (int a = 0; a < 9; a++) {
			//	d2r_dx2.push_back(d2rvecr_dr2[a]);
			//}

			//calculation of hessian elements


            for (int k = 0; k < 3; k++) {
			    for (int l = 0; l < 3; l++) {

					const double hessianValue = dE_dr * d2rvecr_dr2[3 * k + l] + d2E_dr2 * dvecr_dr[k] * dvecr_dr[l];
					cout << "term 1 =" << dE_dr << ", " << d2rvecr_dr2[3 * k + l] << endl;
					cout << "term 2 =" << d2E_dr2 << ", " << dvecr_dr[k] << ", " << dvecr_dr[l] << endl;
					cout << "hessianValue = " << hessianValue << endl;
					
					hessianMatrix[3 * i + k][3 * i + l] += hessianValue;
					hessianMatrix[3 * i + k][3 * j + l] -= hessianValue;
					hessianMatrix[3 * j + k][3 * i + l] -= hessianValue;
					hessianMatrix[3 * j + k][3 * j + l] += hessianValue;

				    //hessianMatrix[3 * i + k][3 * i + l] += dE_dr * d2r_dx2[3 * k + l] + d2E_dr2 * dr_dx[k] * dr_dx[l]; cout << dE_dr * d2r_dx2[3 * k + l] + d2E_dr2 * dr_dx[k] * dr_dx[l] << " " << 3 * i + k << " " << 3 * i + l << endl;
					////hessianMatrix[3 * j + k][3 * j + l] += dE_dr * d2r_dx2[3 * k + l] + d2E_dr2 * dr_dx[3 + k] * dr_dx[3 + l];
					//hessianMatrix[3 * j + k][3 * j + l] += 999999;
					////hessianMatrix[3 * i + k][3 * j + l] += dE_dr * d2r_dx2[9 + (3 * k + l)] + d2E_dr2 * dr_dx[k] * dr_dx[3 + l];
					//hessianMatrix[3 * i + k][3 * j + l] += 999999;
					////hessianMatrix[3 * j + k][3 * i + l] += dE_dr * d2r_dx2[9 + (3 * k + l)] + d2E_dr2 * dr_dx[3 + k] * dr_dx[l];
					//hessianMatrix[3 * j + k][3 * i + l] += 999999;
				}
			}
		}
	}
	return hessianMatrix;
}

int diag (const vector< vector<double> > &hessian) {

	const int hessianSize = hessian.size();
	double hessianArray[hessianSize * hessianSize];

	for (int i = 0; i < hessianSize; i++) {
		for (int j = 0; j < hessianSize; j++) {
			hessianArray[3 * i + j] = hessian[i][j];
		}
	}

	gsl_matrix_view m = gsl_matrix_view_array (hessianArray, hessianSize, hessianSize);

	gsl_vector *eval = gsl_vector_alloc (hessianSize);

	gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc (hessianSize);

	gsl_eigen_symm (&m.matrix, eval, w);

	gsl_eigen_symm_free (w);

    for (int i = 0; i < hessianSize; i++) {
		double eval_i = gsl_vector_get (eval, i);

		cout << eval_i << endl;
	}

	gsl_vector_free (eval);
	return 0;
}

double structure::sumOverAllInteractions (const vector<double> &p) {
    double totalEnergy = 0;
    //iterate over double index ij, where N>j>i and N>i
	for (structure::const_iterator iter = this->begin(); iter != this->end(); ++iter) { 
		for (structure::const_iterator jter = iter + 1; jter != this->end(); ++jter) {
			//cout << "iter is " << *iter << endl;
			//cout << "jter is " << *jter << endl;
            totalEnergy += LJEnergy (coord3d::dist (*iter,*jter), p[0], p[1]);
        }
	}
	return totalEnergy;
}


vector<coord3d> structure::sumOverAllGradients (const vector<double> &p) {
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
	vector<double> *p = static_cast<vector<double>* >(params);
	for (structure::const_iterator iter = kissingSphere.begin(); iter != kissingSphere.end(); ++iter) { 
		for (structure::const_iterator jter = iter + 1; jter != kissingSphere.end(); ++jter) {
			//cout << "iter is " << *iter << endl;
			//cout << "jter is " << *jter << endl;
            totalEnergy += LJEnergy (coord3d::dist (*iter,*jter), (*p)[0], (*p)[1]);
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
	vector<double> *p = static_cast<vector<double>*>(params);
	vector<coord3d> gradients = kissingSphere.sumOverAllGradients(*p);
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
structure structure::optimize (const int &algo_switch, const int &potential_switch, vector<double> parameters, const vector<double> opt, vector<double> &allEnergies) {
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
			cerr << "Error, probably bad input of pot" << endl;
			newGeometry.empty();
			return newGeometry;
    }

    min_function.n = (this->size()) * 3;
	min_function.params = static_cast<void*>(&parameters);

	x = gsl_vector_alloc ( (this->size()) *3 );
    for (structure::size_type i = 0; i < this->size(); ++i) {
	    for (int j = 0; j<=2; ++j) {
	        gsl_vector_set(x, 3 * i + j, (*this)[i][j]);
		}
	}


	switch (algo_switch) {
		case 1:
	        T = gsl_multimin_fdfminimizer_vector_bfgs2;
			break;
		default:
			cerr << "Error, probably bad input of opt algo" << endl;
			newGeometry.empty();
			return newGeometry;
	}

	s = gsl_multimin_fdfminimizer_alloc (T, (this->size()) * 3);
    
	double stepSize = opt[2];
	double accuracy = opt[0];
	gsl_multimin_fdfminimizer_set (s, &min_function, x, stepSize, accuracy);

	size_t i = 0;
	do {
		i++;
		status = gsl_multimin_fdfminimizer_iterate (s);
		
		if (status) {
			cerr << "Something went wrong!" << endl;
			cerr << status << endl;
			break;
		}
        
		double absoluteTolerance = opt[1];
		status = gsl_multimin_test_gradient (s->gradient, absoluteTolerance);

		if (status == GSL_SUCCESS) {
			cout << "Minimum found at:\n" << endl;

		    //create structure for optimized geometry
    	    for (size_t i = 0; i < x->size / 3; ++i) {
                coord3d sphere(gsl_vector_get (s->x, 3 * i), gsl_vector_get (s->x, 3 * i + 1), gsl_vector_get (s->x, 3 * i + 2));
			    cout << sphere <<  endl;
    		    newGeometry.push_back(sphere);
    	    }
        allEnergies.push_back(s->f);
		}
	}
	while (status == GSL_CONTINUE && i <= nsteps);

    
    cout << "-----------------------------------------------" << endl;
    cout << "LJ Energy is: " << s->f << endl;
    if (status == GSL_SUCCESS) {
		cout << "Optimization successful!" << endl;
	}
	else {
		cout << "Error, " << status << ". Optimization failed." << endl;
	}
    gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);

    return newGeometry;
}
