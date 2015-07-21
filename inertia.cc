#include "spheres.h"

using namespace std;

coord3d structure::centreOfMass () {
	double totalMass = nAtoms(); 
	coord3d coordinateSum;
	for (int i = 0; i < nAtoms(); i++) {
		coordinateSum += (*this)[i];
	}
	return coordinateSum * (1 / totalMass);
}


void structure::shiftToCoM (coord3d &CoM) {
	vector<coord3d> shiftedCoordinates;
	for (int i = 0; i < nAtoms(); i++) {
		shiftedCoordinates.push_back((*this)[i] - CoM);
	}
	setCoordinates(shiftedCoordinates);
}

vector< vector<double> > structure::momentOfInertia () {
	vector< vector<double> > inertiaMatrix (3, vector<double> (3, 0));
	for (int i = 0; i < nAtoms(); i++) {
		inertiaMatrix[0][0] += pow((*this)[i][1],2) + pow((*this)[i][2],2);
		inertiaMatrix[1][1] += pow((*this)[i][2],2) + pow((*this)[i][0],2);
		inertiaMatrix[2][2] += pow((*this)[i][0],2) + pow((*this)[i][1],2);
		inertiaMatrix[0][1] += - (*this)[i][0] * (*this)[i][1];
		inertiaMatrix[0][2] += - (*this)[i][0] * (*this)[i][2];
		inertiaMatrix[1][2] += - (*this)[i][1] * (*this)[i][2];
		inertiaMatrix[1][0] += - (*this)[i][0] * (*this)[i][1];
		inertiaMatrix[2][0] += - (*this)[i][0] * (*this)[i][2];
		inertiaMatrix[2][1] += - (*this)[i][1] * (*this)[i][2];
	}
	return inertiaMatrix;
}
		
