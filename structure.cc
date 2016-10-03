#include "structure.h"
#include "lina.h"

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

		
matrix3d structure::m3d_momentOfInertia () {
	matrix3d inertiaMatrix;
	for (int i = 0; i < nAtoms(); i++) {
		inertiaMatrix(0,0) += pow((*this)[i][1],2) + pow((*this)[i][2],2);
		inertiaMatrix(1,1) += pow((*this)[i][2],2) + pow((*this)[i][0],2);
		inertiaMatrix(2,2) += pow((*this)[i][0],2) + pow((*this)[i][1],2);
		inertiaMatrix(0,1) += - (*this)[i][0] * (*this)[i][1];
		inertiaMatrix(0,2) += - (*this)[i][0] * (*this)[i][2];
		inertiaMatrix(1,2) += - (*this)[i][1] * (*this)[i][2];
		inertiaMatrix(1,0) += - (*this)[i][0] * (*this)[i][1];
		inertiaMatrix(2,0) += - (*this)[i][0] * (*this)[i][2];
		inertiaMatrix(2,1) += - (*this)[i][1] * (*this)[i][2];
	}
	return inertiaMatrix;
}

matrix3d structure::m3d_principalAxis () {
	matrix3d I = this->m3d_momentOfInertia();
	matrix3d V = m3d_diagv (I);
	return V;
}

void structure::rotateToPrincipalAxis (matrix3d &principalAxis) {
	matrix3d transMatrix = principalAxis.inverse();
	//matrix3d test =  principalAxis * transMatrix;
	vector<coord3d> coord = this->getCoordinates(), newCoord;

	for (vector<coord3d>::size_type i = 0; i < coord.size(); i++) {
		newCoord.push_back (transMatrix * coord[i]);
	}

	structureCoordinates = newCoord;	
}	

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


bool structure::isMinimum () {
	vector<double> hessian = this->getHessian();
	if (hessian[0] > -0.001) return true;
	return false;
}

vector< vector<int> > structure::createAdjMatrix (vector<double> &p) {
	vector< vector<int> > adjMatrix;
	vector<coord3d> currentCoord = this->getCoordinates();
	for (vector<coord3d>::size_type j = 0; j < currentCoord.size(); j++) {
		vector<int> currentRow;
		for (vector<coord3d>::size_type k = 0; k < currentCoord.size(); k++) {
			if (j == k) currentRow.push_back(1);
			else {
				double diff = fabs(coord3d::dist (currentCoord[j], currentCoord[k]) - p[1]);
				if (diff < 0.1) currentRow.push_back(1);
				else currentRow.push_back(0);
			}
		}
		adjMatrix.push_back(currentRow);
	}
	return adjMatrix;
}

vector<int> structure::createBondVector () {
	vector< vector<int> > adjMatrix = this->getAdjMatrix();
	vector<int> bondVector (adjMatrix.size());
	for (vector< vector<int> >::size_type i = 0; i < adjMatrix.size(); i++) {
		for (vector< vector<int> >::size_type j = 0; j < adjMatrix.size(); j++) {
			bondVector[i] += adjMatrix[i][j] * adjMatrix[j][i];
		}
	}
	return bondVector;
}


vector<double> structure::createAdjMatrix_egenvalues () {
	vector< vector<int> > adjMatrix = this->getAdjMatrix();
	vector< vector<double> > adjMatrix_double;
	for (vector< vector<int> >::size_type i = 0; i < adjMatrix.size(); i++) {
		vector<double> row (adjMatrix[i].begin(), adjMatrix[i].end());
		adjMatrix_double.push_back(row);
	}
	vector<double> adjMatrix_eigenvalues = diag (adjMatrix_double);
	return adjMatrix_eigenvalues;
}

bool structure::compareCoordinates (structure &y) const {
	vector<coord3d> a = this->getCoordinates();
	vector<coord3d> b = y.getCoordinates();

	for (vector<coord3d>::size_type i = 0; i < a.size(); i++) {
		bool match(0);
		for (vector<coord3d>::size_type j = 0; j < b.size(); j++) {
			coord3d diff = a[i] - b[j];
			//cout << i << " " << diff.norm() << endl;
			if (diff.norm() < 1e-4) {
				match = 1;
				break;
			}
		}
		if (match == 0) return false;
	}
	return true;
}

//symmetry
//sig mirror; i: 0=yz 1=xz, 2=xy
vector<coord3d> structure::sig (int a) {
	vector<coord3d> coord = structureCoordinates;
	for (int i = 0; i < this->nAtoms(); i++) { 
		coord[i][a] = - structureCoordinates[i][a];
	}
	return coord;
}

//c2 rotation; ab: 12=x, 02=y, 01=z
vector<coord3d> structure::c2 (int a, int b) {
	vector<coord3d> coord = structureCoordinates;
	for (int i = 0; i < this->nAtoms(); i++) { 
		coord[i][a] = -structureCoordinates[i][a];
		coord[i][b] = -structureCoordinates[i][b];
	}
	return coord;
}

//inversion
vector<coord3d> structure::inv () {
	vector<coord3d> coord = structureCoordinates;
	for (int i = 0; i < this->nAtoms(); i++) { 
		coord[i] = - structureCoordinates[i];
	}
	return coord;
}
