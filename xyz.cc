#include <iostream>
#include <fstream>
#include <iomanip>
#include "spheres.h"

using namespace std;

#define container_output(container) \
template <typename T> ostream& operator<<(ostream& s, const container<T>& v) \
	{ \
	s << "{"; \
	for(typename container<T>::const_iterator x(v.begin());x!=v.end();){ \
		s << *x; \
		if(++x!=v.end()) s << ","; \
	} \
	s << "}"; \
	return s; \
	}
container_output(vector);


void xyzout (structure &outputStructure, const string &name = "structure.xyz") {
	vector<coord3d> coordinates = outputStructure.getCoordinates();
	ofstream xyz;
	xyz.open (name);
	xyz << outputStructure.nAtoms() << endl;
	xyz << "Energy: " << outputStructure.getEnergy() << ", moment of inertia: " << outputStructure.getMomentOfInertia() << endl;
	for (int i = 0; i < outputStructure.nAtoms(); i++) {
		xyz << setw(5) << left << "X" << setw(15) << setprecision(6) << setiosflags(ios::fixed) << right << outputStructure[i][0] << setw(15) << outputStructure[i][1] << setw(15) << outputStructure[i][2] << endl;
	}
	xyz.close();
}
	
