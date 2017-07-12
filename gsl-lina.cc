#include <iostream>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <algorithm>
#include <assert.h>
#include "geometry.h"


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
        cout << "An error occured, status " << status << ", " << gsl_strerror(status) << endl;
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

vector<pair<double, vector<double> > > diagv (vector< vector<double> > &matrix) {

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
        cout << "An error occured, status " << status << ", " << gsl_strerror(status) << endl;
    }


    gsl_eigen_symmv_free (w);
    
    vector<pair<double, vector<double> > > eval_evec;
    for (int i = 0; i < matrixSize; i++) {
        gsl_vector_view evec_i = gsl_matrix_column (evec, i);
        vector<double> eigv;
        for (int j = 0; j < matrixSize; j++) {
            eigv.push_back(gsl_vector_get (&evec_i.vector, j));
        }
        eval_evec.push_back(make_pair(gsl_vector_get (eval, i), eigv));
    }


    gsl_vector_free (eval);
    gsl_matrix_free (evec);
    return eval_evec;
}

matrix3d m3d_diagv (matrix3d &matrix) {


    gsl_matrix_view m = gsl_matrix_view_array (matrix.values, 3, 3);
    gsl_vector *eval = gsl_vector_alloc (3);
    gsl_matrix *evec = gsl_matrix_alloc (3, 3);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (3);

    int status = gsl_eigen_symmv (&m.matrix, eval, evec, w);

    if (status) {
        cout << "An error occured, status " << status << ", " << gsl_strerror(status) << endl;
    }

    gsl_eigen_symmv_free (w);


    vector<pair<double, vector<double> > > eval_evec;
    for (int i = 0; i < 3; i++) {
        gsl_vector_view evec_i = gsl_matrix_column (evec, i);
        vector<double> eigv;
        for (int j = 0; j < 3; j++) {
            eigv.push_back(gsl_vector_get (&evec_i.vector, j));
        }
        eval_evec.push_back(make_pair(gsl_vector_get (eval, i), eigv));
    }
    
    auto pairCompare = [&] (const pair<double, vector<double> > a, const pair<double, vector<double> > b) {
                return a.first < b.first;
    };
    sort (eval_evec.begin(), eval_evec.end(), pairCompare); 

    matrix3d diag;

    for (int i = 0; i < 3; i++) {
        diag(0,i) = eval_evec[i].second[0]; 
        diag(1,i) = eval_evec[i].second[1]; 
        diag(2,i) = eval_evec[i].second[2];
    }

    gsl_vector_free (eval);
    gsl_matrix_free (evec);

    return diag;
}
