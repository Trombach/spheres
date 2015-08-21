#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <algorithm>
#include <assert.h>
#include "spheres.h"


using namespace std;

vector< vector<double> > structure::hessian (const vector<double> &p) {
	const double epsilon = p[0];
	const double rm = p[1];
	const double exp1 = p[2];
	const double exp2 = p[3];
	vector < vector<double>  > hessianMatrix (this->nAtoms() * 3, vector<double> (this->nAtoms() * 3, 0));
    for (int i = 0; i < this->nAtoms(); i++) {
		for (int j = i + 1; j < this->nAtoms(); j++) {

			const coord3d vecr = (*this)[i] - (*this)[j];
			const double r = coord3d::dist((*this)[i], (*this)[j]);

			//calculate first derivative
			double dE_dr = - ( epsilon / rm ) * ( exp1 * (pow (rm / r, exp1 + 1)) - 2 * exp2 * (pow (rm / r, exp2 + 1)) );

			//calculate second derivative
            double d2E_dr2 = epsilon / pow (rm, 2) * ( (pow (exp1, 2) + exp1) * pow (rm / r, exp1 + 2) - (2 * pow (exp2, 2) + 2 * exp2) * pow (rm / r, exp2 + 2) );

			//calculate derivatives of r
			coord3d dvecr_dr = coord3d::dnorm(vecr);
			//cout << "dvecr_dr" << endl;
			//for (int i = 0; i < 3; i++) {cout << dvecr_dr[i] << " ";}
			//cout << endl << "------------------" <<endl;

			vector<double> d2rvecr_dr2(9, double());
			coord3d::ddnorm(vecr, d2rvecr_dr2);


			//calculation of hessian elements

			//loop over all 6 coordinates of 1 atom pair
            for (int k = 0; k < 3; k++) {
			    for (int l = 0; l < 3; l++) {
					//calculate the value first, which will always only differ by sign
					const double hessianValue = dE_dr * d2rvecr_dr2[3 * k + l] + d2E_dr2 * dvecr_dr[k] * dvecr_dr[l];
					//cout << "term 1 =" << dE_dr << ", " << d2rvecr_dr2[3 * k + l] << endl;
					//cout << "term 2 =" << d2E_dr2 << ", " << dvecr_dr[k] << ", " << dvecr_dr[l] << endl;
					//cout << "hessianValue = " << hessianValue << endl;
					
					//write hessian
					//this is basically a 2 atom hessian, where the diagonal quadrants are the same and the remaining quadrants are of the opposite sign
					//each quadrant can have contributions from different atom pairs
					//eg atom pair 1/2 and 1/3 both have non zero second derivaties with respect to the coordinates of atom 1
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

vector<double> diag (vector< vector<double> > &matrix) {

	for (vector< vector<double> >::size_type i = 0; i < matrix.size(); i++) {
		assert (matrix.size() == matrix[i].size());
	}
	
	const int matrixSize = matrix.size();
	double matrixArray[matrixSize * matrixSize];

	for (int i = 0; i < matrixSize; i++) {
		for (int j = 0; j < matrixSize; j++) {
			matrixArray[matrixSize * i + j] = matrix[i][j];
		}
	}
	

	gsl_matrix_view m = gsl_matrix_view_array (matrixArray, matrixSize, matrixSize);
	gsl_vector *eval = gsl_vector_alloc (matrixSize);
	gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc (matrixSize);

	int status = gsl_eigen_symm (&m.matrix, eval, w);

	if (status) {
		cout << "An error occured, status " << status << endl;
	}

	gsl_eigen_symm_free (w);
    
	vector<double> eigenvalues;
    for (int i = 0; i < matrixSize; i++) {
		eigenvalues.push_back(gsl_vector_get (eval, i));
		sort(eigenvalues.begin(), eigenvalues.end());
	}

	gsl_vector_free (eval);
	return eigenvalues;
}

pair<vector<double>, vector< vector<double> > > diagv (vector< vector<double> > &matrix) {

	for (vector< vector<double> >::size_type i = 0; i < matrix.size(); i++) {
		assert (matrix.size() == matrix[i].size());
	}
	
	const int matrixSize = matrix.size();
	double matrixArray[matrixSize * matrixSize];

	for (int i = 0; i < matrixSize; i++) {
		for (int j = 0; j < matrixSize; j++) {
			matrixArray[matrixSize * i + j] = matrix[i][j];
		}
	}
	

	gsl_matrix_view m = gsl_matrix_view_array (matrixArray, matrixSize, matrixSize);
	gsl_vector *eval = gsl_vector_alloc (matrixSize);
	gsl_matrix *evec = gsl_matrix_alloc (matrixSize, matrixSize);
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (matrixSize);

	int status = gsl_eigen_symmv (&m.matrix, eval, evec, w);

	if (status) {
		cout << "An error occured, status " << status << endl;
	}


	gsl_eigen_symmv_free (w);
    
	vector<double> eigenvalues;
    for (int i = 0; i < matrixSize; i++) {
		eigenvalues.push_back(gsl_vector_get (eval, i));
	}
	vector< vector<double> > eigenvectors;
	for (int i = 0; i < matrixSize; i++) {
		vector<double> eigv;
		for (int j = 0; j < matrixSize; j++) {
			eigv.push_back(gsl_matrix_get(evec, i, j));
		}
		eigenvectors.push_back(eigv);
	}

	gsl_vector_free (eval);
	gsl_matrix_free (evec);
	return make_pair(eigenvalues,eigenvectors);
}
