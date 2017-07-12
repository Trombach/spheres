#include <iostream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <dlib/matrix.h>
#include "geometry.h"

#define DLIB_USE_LAPACK

using namespace std;

typedef dlib::matrix<double> dlib_matrix;
typedef dlib::matrix<double,3,3> dlib_matrix3d;


vector<double> diag (vector< vector<double> > &matrix) 
{

    for (vector< vector<double> >::size_type i = 0; i < matrix.size(); i++) 
    {
        assert (matrix.size() == matrix[i].size());
    }
    
    const int matrixSize = matrix.size();
    dlib_matrix matrixArray;
    matrixArray.set_size(matrixSize, matrixSize);

    for (int i = 0; i < matrixSize; i++) 
    {
        for (int j = 0; j < matrixSize; j++) 
        {
            matrixArray(i, j) = matrix[i][j];
        }
    }
    
    dlib::eigenvalue_decomposition<dlib_matrix> eigen (matrixArray);
    dlib::eigenvalue_decomposition<dlib_matrix>::column_vector_type eval = eigen.get_real_eigenvalues();

    
    vector<double> eigenvalues;
    for (int i = 0; i < matrixSize; i++) 
    {
        eigenvalues.push_back(eval(i));
    }
    sort(eigenvalues.begin(), eigenvalues.end());

    return eigenvalues;
}

vector<pair<double, vector<double> > > diagv (vector< vector<double> > &matrix) 
{

    for (vector< vector<double> >::size_type i = 0; i < matrix.size(); i++) 
    {
        assert (matrix.size() == matrix[i].size());
    }
    
    const int matrixSize = matrix.size();
    dlib_matrix matrixArray;
    matrixArray.set_size(matrixSize, matrixSize);

    for (int i = 0; i < matrixSize; i++) 
    {
        for (int j = 0; j < matrixSize; j++) 
        {
            matrixArray(i, j) = matrix[i][j];
        }
    }
    
    dlib::eigenvalue_decomposition<dlib_matrix> eigen (matrixArray);
    dlib::eigenvalue_decomposition<dlib_matrix>::column_vector_type eval = eigen.get_real_eigenvalues();    
    dlib::eigenvalue_decomposition<dlib_matrix>::matrix_type evec = eigen.get_pseudo_v();

    
    vector<pair<double, vector<double> > > eval_evec;
    for (int i = 0; i < matrixSize; i++) 
    {
        vector<double> eigv;
        for (int j = 0; j < matrixSize; j++) 
        {
            eigv.push_back(evec(j, i));
        }
        eval_evec.push_back(make_pair(eval(i), eigv));
    }


    return eval_evec;
}

matrix3d m3d_diagv (matrix3d &matrix) {

    dlib_matrix3d matrixArray;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            matrixArray(i,j) = matrix(i,j);
        }
    }

    dlib::eigenvalue_decomposition<dlib_matrix3d> eigen (matrixArray);
    dlib::eigenvalue_decomposition<dlib_matrix3d>::column_vector_type eval = eigen.get_real_eigenvalues();
    dlib::eigenvalue_decomposition<dlib_matrix3d>::matrix_type evec = eigen.get_pseudo_v();
    

    vector<pair<double, vector<double> > > eval_evec;
    for (int i = 0; i < 3; i++) {
        vector<double> eigv;
        for (int j = 0; j < 3; j++) {
            eigv.push_back(evec(j, i));
        }
        eval_evec.push_back(make_pair(eval(i), eigv));
    }
    
    auto pairCompare = [&] (const pair<double, vector<double> > a, const pair<double, vector<double> > b) 
    {
                return a.first < b.first;
    };
    sort (eval_evec.begin(), eval_evec.end(), pairCompare); 

    matrix3d diag;

    for (int i = 0; i < 3; i++) {
        diag(0,i) = eval_evec[i].second[0]; 
        diag(1,i) = eval_evec[i].second[1]; 
        diag(2,i) = eval_evec[i].second[2];
    }


    return diag;
}
