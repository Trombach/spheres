#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include "iop.h"

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
	xyz.open ("output/" + name);
	xyz << outputStructure.nAtoms() << endl;
	xyz << "Energy: " << outputStructure.getEnergy() << ", moment of inertia: " << outputStructure.getMomentOfInertia() << endl;
	for (int i = 0; i < outputStructure.nAtoms(); i++) {
		xyz << setw(5) << left << "X" << setw(15) << setprecision(6) << setiosflags(ios::fixed) << right << outputStructure[i][0] << setw(15) << outputStructure[i][1] << setw(15) << outputStructure[i][2] << endl;
	}
	xyz.close();
}


void simpleout (vector<structure> &outputStructures, stringstream &output) {
	for (vector<structure>::size_type i=0; i < outputStructures.size(); i++) {
		for (int j=0; j < outputStructures[i].nAtoms(); j++) {
			output << outputStructures[i][j][0] << " " << outputStructures[i][j][1] << " " << outputStructures[i][j][2] << endl;
		}
		output << endl;
	}
}



//function for determining if a line in the input is emtpy
//used for breaking the read loop
bool justempty(string str) {
	if (str.empty()) {
		return true;
	}
	return false;
}

//function for testing if input file exists
//returns true for existing file
bool fexists (const std::string &fileName) {
    ifstream infile(fileName.c_str());
	bool exist = infile.good();
	infile.close();
    return exist;
}

//function for reading all strucutures of input file
//structures must be separated by a blank line, last line should be a blank line as well
vector<structure> readallstruct (const std::string& fileName) {
	ifstream infile(fileName.c_str());
	coord3d sphere; //sphere refers to 1 'atom'
	string line;
	int number = 1;
	vector<structure> allKissingSpheres; //this vector stores all read in structures
	//double while loop, inner loop goes over all coordinates of one structure, outer loop goes over all 
	//structures till eof
	while (true) {
	    structure kissingSphere; //kissingSphere refers to 1 cluster of kissing Spheres
		vector<coord3d> coordinates;
	    while (true) {
	        getline(infile,line);
	        if(justempty(line) || infile.eof()) break;
	        else {
	    	    stringstream lineStream(line); //string can't be read into coord3d object, convert to 
				                               //stringstream first
	    	    lineStream >> sphere;
				coordinates.push_back(sphere);
	        }
	    }
		kissingSphere.setCoordinates(coordinates);
		kissingSphere.setNumber(number);
		if (infile.eof()) break; //does last line has to be blank? break works for files with blank line at 
		                         //the end
		else {
		    allKissingSpheres.push_back(kissingSphere);
			number++;
	    }
	}
    infile.close();
	cout << "\t" << "Number of structures: " << allKissingSpheres.size() << endl;
	return allKissingSpheres;
}
