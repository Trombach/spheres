#include "geometry.h"


matrix3d matrix3d::transpose() const {
    const matrix3d &M(*this);
    matrix3d Mt;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            Mt(i,j) = M(j,i);
    return Mt;
}

double matrix3d::norm() const {
    const matrix3d &M(*this);

    double norm = 0;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            norm += M(i,j)*M(i,j);

    return sqrt(norm);
}

double matrix3d::det() const {
    const matrix3d &M(*this);

    double det = 0;

    det += M(0,0) * M(1,1) * M(2,2);
    det += M(0,1) * M(1,2) * M(2,0);
    det += M(0,2) * M(1,0) * M(2,1);

    det -= M(0,2) * M(1,1) * M(2,0);
    det -= M(0,1) * M(1,0) * M(2,2);
    det -= M(0,0) * M(1,2) * M(2,1);

    return det;
}

matrix3d matrix3d::inverse() {
    matrix3d &M (*this);
    matrix3d Mi;

    double det = M.det();   

    Mi(0,0) = M(1,1) * M(2,2) - M(1,2) * M(2,1);
    Mi(0,1) = - (M(1,0) * M(2,2) - M(1,2) * M(2,0));
    Mi(0,2) = M(1,0) * M(2,1) - M(1,1) * M(2,0);
    Mi(1,0) = - (M(0,1) * M(2,2) - M(0,2) * M(2,1));
    Mi(1,1) = M(0,0) * M(2,2) - M(0,2) * M(2,0);
    Mi(1,2) = - (M(0,0) * M(2,1) - M(0,1) * M(2,0));
    Mi(2,0) = M(0,1) * M(1,2) - M(0,2) * M(1,1);
    Mi(2,1) = - (M(0,0) * M(1,2) - M(0,2) * M(1,0));
    Mi(2,2) = M(0,0) * M(1,1) - M(0,1) * M(1,0);

    return Mi.transpose() * (1/det);
}


matrix3d matrix3d::operator*(const matrix3d& B) const {
    const matrix3d &A(*this);
    matrix3d C;

    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            double sum = 0;
            for(int k=0;k<3;k++){
                sum += A(i,k)*B(k,j);
            }
            C(i,j) = sum;
        }
    }
    return C;
}


coord3d matrix3d::operator*(const coord3d& x) const {
    coord3d y;
    for(int j=0;j<3;j++)
        y += coord3d(values[j]*x[j],values[3+j]*x[j],values[6+j]*x[j]);
    return y;
}


std::vector<coord3d> matrix3d::operator*(const std::vector<coord3d>& xs) const {
    const matrix3d &A(*this);
    std::vector<coord3d> ys(xs.size());

    for(int i=0;i<xs.size();i++) ys[i] = A*xs[i];
    return ys;
}


//structure matrix3d::operator*(structure &xs) const {
//const matrix3d &A(*this);
//structure ys;

//for (int i = 0; i < xs.nAtoms(); i++) ys.push_back(A * xs[i]);
//return ys;
//} 
