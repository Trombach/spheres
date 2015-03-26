#ifndef SPHERES
#define SPHERES

#include <iostream>
#include <vector>
#include <math.h>


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
    

public:   

    //
    //function to sum over all sphere interactions, change later to work with different potentials
    //
   	double sumOverAllInteractions ();

	//
	//function to sum over all gradients to get gradient for each sphere
	//
	vector<coord3d> sumOverAllGradients (const vector<double> &p);
    
	//
	//initialize gsl minimizer function
	//
	int optimize (const int &algo_switch, const int &potential_switch, const vector<double> parameters, const vector<double> opt);
};

#endif
