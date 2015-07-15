#ifndef SPHERES
#define SPHERES

#include <iostream>
#include <vector>
#include <math.h>
#include <string>


std::vector<double> diag (std::vector< std::vector<double> > &hessian); 

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

    friend std::ostream& operator<<(std::ostream &s, const coord3d& x){ s << std::fixed << "{" << x[0] << "," << x[1] << "," << x[2]<< "}"; return s; }
    friend std::istream& operator>>(std::istream &s, coord3d& x){ for(int i=0;i<3;i++){ s >> x[i]; } return s; }

};


class structure:public std::vector<coord3d> {    
    
private:
	double structureEnergy;
	int structureNumber;
	vector<coord3d> structureCoordinates;

public:   
	
	structure() {};
	structure(int number, double energy, vector<coord3d> coordinates) { structureNumber = number; structureEnergy = energy; structureCoordinates = coordinates; }

	int getNumber() const { return structureNumber; }
	double getEnergy() const { return structureEnergy; }
	vector<coord3d> getCoordinates() const { return structureCoordinates; }

	void setNumber (int number) { structureNumber = number; }
	void setEnergy (double energy) { structureEnergy = energy; }
	void setCoordinates (vector<coord3d> coordinates) { structureCoordinates = coordinates; }

	structure &operator*= (const double &y) {
		for (vector<coord3d>::size_type i = 0; i < this->size(); i++) { (*this)[i] *= y; }
		return *this;
	}
	structure operator* (const double &y) const { return structure(*this) *= y; }

    //function to sum over all sphere interactions, change later to work with different potentials
   	double sumOverAllInteractions (const vector<double> &p);
	//function to sum over all gradients to get gradient for each sphere
	vector<coord3d> sumOverAllGradients (const vector<double> &p);
	//initialize gsl minimizer function
	structure optimize (const int &algo_switch, const int &potential_switch, const vector<double> parameters, const vector<double> opt, vector<double> &allEnergies);

	vector< vector<double> > hessian (const vector<double> &p);

};
    
#endif
