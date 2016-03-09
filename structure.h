#ifndef STRUCTURE
#define STRUCTURE

#include "geometry.h"

class structure {    
    
private:
	double structureEnergy;
	int structureNumber;
	std::vector<coord3d> structureCoordinates;
	std::vector<double> structureMomentOfInertia;
	std::vector<double> structureHessian;

public:   
	structure() {};
	structure(int number, double energy, std::vector<coord3d> coordinates) { structureNumber = number; structureEnergy = energy; structureCoordinates = coordinates; }

	int getNumber() const { return structureNumber; }
	double getEnergy() const { return structureEnergy; }
	std::vector<coord3d> getCoordinates() const { return structureCoordinates; }
	std::vector<double> getMomentOfInertia() const { return structureMomentOfInertia; }
	std::vector<double> getHessian() const { return structureHessian; }

	void setNumber (int number) { structureNumber = number; }
	void setEnergy (double energy) { structureEnergy = energy; }
	void setCoordinates (std::vector<coord3d> coordinates) { structureCoordinates = coordinates; }
	void setMomentOfInertia (std::vector<double> inertiaEigenvalues) { structureMomentOfInertia = inertiaEigenvalues; }
	void setHessian (std::vector<double> hessianEigenvalues) {structureHessian = hessianEigenvalues; }

	int nAtoms() { return (this->getCoordinates()).size(); }

	void push_back(coord3d spheres) { structureCoordinates.push_back(spheres); }

	structure &operator*= (const double &y) {
		for (int i = 0; i < this->nAtoms(); i++) { (*this)[i] *= y; }
		return *this;
	}
	structure operator* (const double &y) const { return structure(*this) *= y; }
	coord3d& operator[] (unsigned int i) { return structureCoordinates[i]; }
	bool operator< (const structure y) const { return this->getEnergy() < y.getEnergy(); }

    //function to sum over all sphere interactions, change later to work with different potentials
   	double sumOverAllInteractions (const std::vector<double> &p);
	//function to sum over all gradients to get gradient for each sphere
	std::vector<coord3d> sumOverAllGradients (const std::vector<double> &p);
	//initialize gsl minimizer function
	structure optimize (std::ofstream &min, const int &algo_switch, const int &potential_switch, const std::vector<double> parameters, const std::vector<double> opt, std::vector<double> &allEnergies);

	std::vector< std::vector<double> > hessian (const std::vector<double> &p);

	coord3d centreOfMass ();
	void shiftToCoM (coord3d &CoM);
	std::vector< std::vector<double> > momentOfInertia ();
	bool isMinimum ();


};

#endif
