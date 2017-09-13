#ifndef GEOMETRY
#define GEOMETRY

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <memory.h>


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
    // d/dx_i ||x|| = x_i/||x||.
    static coord3d dnorm(const coord3d& x){ return x/x.norm(); }
    // d^2/(dx_i dx_j) ||x|| = -x_i x_j/||x||^3 + [i==j]/||x||
    static void ddnorm(const coord3d& x, std::vector<double> &H)
    {
      const double n = 1.0/x.norm(), n3 = n*n*n;

      for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            H[i*3+j] = -x[i]*x[j]*n3 + (i==j? n : 0);
    }

    friend std::ostream& operator<<(std::ostream &s, const coord3d& x){ s << std::fixed << x[0] << " " << x[1] << " " << x[2]; return s; }
    friend std::istream& operator>>(std::istream &s, coord3d& x){ for(int i=0;i<3;i++){ s >> x[i]; } return s; }

};


struct matrix3d {
    double values[9];

//  matrix3d()                { memset(values,0,9*sizeof(double)); }
    matrix3d(const double *v) { memcpy(values,v,9*sizeof(double)); }
    explicit matrix3d(const double r=0, const double s=0, const double t=0, const double u=0, const double v=0, const double w=0, const double x=0, const double y=0, const double z=0) {
        values[0]=r; values[1]=s; values[2]=t; values[3]=u; values[4]=v; values[5]=w; values[6]=x; values[7]=y; values[8]=z;
    }

    double& operator()(int i, int j)       { return values[i*3+j]; }
    double  operator()(int i, int j) const { return values[i*3+j]; }
    matrix3d operator+(const matrix3d& y) const { return matrix3d(*this) += y; }
    matrix3d operator-(const matrix3d& y) const { return matrix3d(*this) -= y; }
    matrix3d& operator+=(const matrix3d& y){ for(int i=0;i<3;++i){for(int j=0;j<3;++j){values[3*i+j] += y(i,j);}}; return *this; }
    matrix3d& operator-=(const matrix3d& y){ for(int i=0;i<3;++i){for(int j=0;j<3;++j){values[3*i+j] -= y(i,j);}}; return *this; }
    matrix3d operator*(const double s)   const { return matrix3d(*this) *= s; }
    matrix3d& operator*=(const double& s){ for(int i=0;i<3;++i){for(int j=0;j<3;++j){values[3*i+j] *= s;}}; return *this; }
    matrix3d operator-() const {matrix3d m(-values[0],-values[1],-values[2],-values[3],-values[4],-values[5],-values[6],-values[7],-values[8]); return m;}

    matrix3d transpose() const;

    double norm() const;

    double det() const;

    matrix3d inverse();
  
    matrix3d operator*(const matrix3d& B) const;

    coord3d operator*(const coord3d& x) const;
  
    std::vector<coord3d> operator*(const std::vector<coord3d>& xs) const;
  

    //structure operator*(structure &xs) const;
};  
#endif
