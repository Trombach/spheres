#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <algorithm>
#include <assert.h>


using namespace std;


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
