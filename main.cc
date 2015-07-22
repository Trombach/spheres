#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include "spheres.h"
#include <libconfig.h++>
#include <algorithm>


using namespace std; 
using namespace libconfig;

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

void xyzout (structure &outputStructure, const string &name = "structure.xyz");

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
	cout << "Number of structures is " << allKissingSpheres.size() << endl;
	return allKissingSpheres;
}

//MAIN FUNCTION BEGINS HERE

int main (int argc, char *argv[]) {
	
	//SOME START UP CHECKS AND ARGUMENT PROCESSING
    if ( argc < 2 ) {
        cout << "Please provide filename." << endl;
        return 1;
    }
    else {
        cout << "Filename is " << argv[ argc -1 ] << "."  << endl;
    }

	//CHECK IF FILE EXISTS

    string fileName = argv[ argc -1 ]; //safe input file name
    cout << fileName << endl;
    
    if (fexists(fileName)) {
        cout << "File " << fileName << " exists." << endl;
    }
    else {
        cout << "File " << fileName << " does not exist." << endl;
        return 1;
    }
    if (fexists("settings")) {
		cout << "Settings file found." << endl;
	}
	else {
		cout << "No settings file in working directory" << endl;
		return 1;
	}

	
	//READ SETTINGS FILE

	libconfig::Config cfg;
	try {
	    cfg.readFile("settings");
	}
	catch(const FileIOException &fioex) {
        cout << "I/O error while reading file." << endl;
		return(EXIT_FAILURE);
	}
    catch(const ParseException &pex) {
		cout << "Parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << endl;
	    return(EXIT_FAILURE);
	}

	const Setting &root = cfg.getRoot();

	string potential;
	vector<double> p;
	vector<int> structureNumbers;
	int potential_switch, algo_switch, scaling_switch, output_switch(1);
    string algo;
	double scalingFactor(1.0), accuracy(0.1), dforce(1e-3), stepsize(0.01);
	int nsteps(100);
	
	try {
		const Setting &numbers = root["output"]["number"];
		int count = numbers.getLength();
		for (int i = 0; i < count; i++) {
			structureNumbers.push_back (numbers[i]);
		}
	}
	catch (const SettingNotFoundException &nfex) {
		cout << "Output setting not Found." << endl;
		output_switch = 0;
	}
	cout << "Structures to print: " << structureNumbers << endl;


	if (cfg.lookupValue("scaling.factor", scalingFactor)) {
		cout << "Scaling factor in seetings found. Input coordinates will be scaled by " << scalingFactor << endl;
		scaling_switch = 1;
	}
	else {
		cout << "No scaling of input coordinates." << endl;
	}

    if (cfg.lookupValue("potential.name", potential)) {
		if (potential == "LJ") {
		    potential_switch = 1;
		}
		cout << "Chosen potential is " << potential << endl;
	}
	else {
		cout << "No 'potential' setting in configuration file." << endl;
		return 1;
	}

	
    if (potential == "LJ") {
		double exp1(12), exp2(6), rm, epsilon;
		cout << "Setting LJ Parameters." << endl;
		if (cfg.lookupValue("potential.epsilon", epsilon)) {
			p.push_back(epsilon);
			cout << "Epsilon set to " << p[0] << endl;
		}
		else {
			cout << "No 'epsilon' in configuration file." << endl;
			return 1;
		}
		if (cfg.lookupValue("potential.rm", rm)) {
			p.push_back(rm);
			cout << "Rm set to " << p[1] << endl;
		}
		else {
			cout << "No 'rm' in configuration file." << endl;
			return 1;
		}
		if (cfg.lookupValue("potential.exp1", exp1)) {
			p.push_back(exp1);
			cout << "First exponent set to " << exp1 << endl;
		}
		else {
			cout << "No first exponent in settings file." << endl;
			cout << "Fall back to " << exp1 << endl;
			p.push_back(exp1);
		}
		if (cfg.lookupValue("potential.exp2", exp2)) {
			p.push_back(exp2);
			cout << "First exponent set to " << exp2 << endl;
		}
		else {
			cout << "No first exponent in settings file." << endl;
			cout << "Fall back to " << exp2 << endl;
			p.push_back(exp2);
		}
	}



	vector<double> opt; //vector of algo settings, 0 == accuracy, 1 == dforce, 2 == stepsize, 3 == nsteps
	if (cfg.lookupValue("opt.name", algo)) {
		if (algo == "BFGS") {
			algo_switch = 1;
		}
		cout << "Optimization algo set to " << algo << endl;
	}
	else {
	    cout << "No 'opt' setting in configuration file." << endl;
		return 1;
	}

	if (cfg.lookupValue("opt.accuracy", accuracy)) {
		cout << "Accuracy set to " << accuracy << endl;
	}
	else {
	    cout << "No 'accuracy' setting in configuration file." << endl;
		cout << "Fall back to " << accuracy << endl;
	}

	if (cfg.lookupValue("opt.dforce", dforce)) {
		cout << "Convergence criterion for gradients set to " << dforce << endl;
	}
	else {
	    cout << "No 'dforce' setting in configuration file." << endl;
	    cout << "Fall back to " << dforce << endl;
	}

	if (cfg.lookupValue("opt.stepsize", stepsize)) {
		cout << "Stepsize set to " << stepsize << endl;
	}
	else {
	    cout << "No 'stepsize' setting in configuration file." << endl;
		cout << "Fall back to " << stepsize << endl;
	}

	if (cfg.lookupValue("opt.nsteps", nsteps)) {
		cout << "Number of steps set to " << nsteps << endl;
	}
	else {
	    cout << "No 'nsteps' setting in configuration file." << endl;
		cout << "Fall back to " << stepsize << endl;
	}
    opt.push_back(accuracy);
	opt.push_back(dforce);
	opt.push_back(stepsize);
	opt.push_back(nsteps);

    
	//READ IN ALL STRUCTURES AT ONCE (AND CALCULATE LJ-ENERGY FOR EACH STRUCTURE)
    vector<structure> allKissingSpheres = readallstruct(fileName);
	
	//if scaling is found in settings file scale all coordinates accordingly
	switch (scaling_switch) {
		case 1:
			for (vector<structure>::size_type i = 0; i < allKissingSpheres.size(); i++) {
				allKissingSpheres[i] *= scalingFactor;
			}
		default:
			break;
	}



    vector<structure> optimizedKissingSpheres;
	vector< vector<double> > hessian;
	vector<double> eigenValues;
	vector<double> allEnergies;

	for (vector<structure>::size_type i = 0; i < allKissingSpheres.size(); i++) {
		cout << "Optimization for structure no " << allKissingSpheres[i].getNumber() << endl;
        optimizedKissingSpheres.push_back( allKissingSpheres[i].optimize( algo_switch, potential_switch, p, opt , allEnergies) );
		optimizedKissingSpheres[i].setNumber( allKissingSpheres[i].getNumber() );
		hessian = optimizedKissingSpheres[i].hessian(p);
		eigenValues = diag(hessian);
		cout << "Eigenvalues of the hessian are:" << endl << eigenValues << endl;
		cout << "###############################################################\n" << endl;
	}

	for (vector<structure>::size_type i = 0; i < optimizedKissingSpheres.size(); i++) {
		coord3d CoM = optimizedKissingSpheres[i].centreOfMass();
		optimizedKissingSpheres[i].shiftToCoM(CoM);
		vector< vector<double> > inertiaTensor = optimizedKissingSpheres[i].momentOfInertia();
		vector<double> inertia = diag(inertiaTensor);
		optimizedKissingSpheres[i].setMomentOfInertia(inertia);
	}


	sort(optimizedKissingSpheres.begin(), optimizedKissingSpheres.end());
    for (vector<structure>::size_type i = 0; i < optimizedKissingSpheres.size(); i++) {
		cout << "Number: " << optimizedKissingSpheres[i].getNumber() << "\t" << "Energy: " << optimizedKissingSpheres[i].getEnergy() << "\t" << "Inertia: " << optimizedKissingSpheres[i].getMomentOfInertia() << endl;
	}

	//OUTPUT CONTROL
	
	switch (output_switch) {
		case 1:
			for (vector<structure>::size_type i = 0; i < optimizedKissingSpheres.size(); i++) {
				for (vector<int>::size_type j = 0; j < structureNumbers.size(); j++) {
					if (optimizedKissingSpheres[i].getNumber() == structureNumbers[j]) {
						xyzout(optimizedKissingSpheres[i], "optStruct" + to_string(optimizedKissingSpheres[i].getNumber()));
					}
				}
			}
		default:
			break;
	}

    return 0; 

    
}
