#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <libconfig.h++>
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

using namespace libconfig;

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
	//cout << "\t" << "Number of structures: " << allKissingSpheres.size() << endl;
	return allKissingSpheres;
}

//function to read settings file
int readsettings (vector<double> &opt, vector<double> &p, vector<int> &switches, double &scalingFactor) {

	libconfig::Config cfg;
	try {
	    cfg.readFile("settings");
	}
	catch(const FileIOException &fioex) {
        cerr << "\tI/O error while reading file." << endl;
		return(EXIT_FAILURE);
	}
    catch(const ParseException &pex) {
		cerr << "\tParse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << endl;
	    return(EXIT_FAILURE);
	}

	const Setting &root = cfg.getRoot();

	string potential, algo;
	double accuracy(0.1), dforce(1e-3), stepsize(0.01);
	int potential_switch, algo_switch, scaling_switch, output_switch(1);
	int nsteps(100);
	//bool printInput(false);
	
	//try {
	//	const Setting &numbers = root["output"]["number"];
	//	int count = numbers.getLength();
	//	for (int i = 0; i < count; i++) {
	//		structureNumbers.push_back (numbers[i]);
	//	}
	//}
	//catch (const SettingNotFoundException &nfex) {
	//	cerr << "\tOutput setting not Found." << endl;
	//	output_switch = 0;
	//}
	//cout << "\tStructures to print: " << structureNumbers;
	//if (cfg.lookupValue ("output.input", printInput)) {
	//	cout << ", including input coordinates" << endl;
	//}
	//else {
	//	cout << ", excluding input coordinates" << endl;
	//}


	cout << endl;

	if (cfg.lookupValue("scaling.factor", scalingFactor)) {
		cout << "\tScaling: " << scalingFactor << endl;
		scaling_switch = 1;
	}
	else {
		cout << "\tNo scaling" << endl;
	}

    if (cfg.lookupValue("potential.name", potential)) {
		if (potential == "LJ") {
		    potential_switch = 1;
		}
		cout << "\tPotential: " << potential << endl;
	}
	else {
		cout << "\tNo 'potential' setting in configuration file." << endl;
		return 1;
	}

	
    if (potential == "LJ") {
		double exp1(12), exp2(6), rm, epsilon;
		if (cfg.lookupValue("potential.epsilon", epsilon)) {
			p.push_back(epsilon);
			cout << "\t\tEpsilon: " << p[0] << endl;
		}
		else {
			cout << "\tNo 'epsilon' in configuration file." << endl;
			return 1;
		}
		if (cfg.lookupValue("potential.rm", rm)) {
			p.push_back(rm);
			cout << "\t\tRm: " << p[1] << endl;
		}
		else {
			cout << "\tNo 'rm' in configuration file." << endl;
			return 1;
		}
		if (cfg.lookupValue("potential.exp1", exp1)) {
			p.push_back(exp1);
			cout << "\t\tExp1: " << exp1 << endl;
		}
		else {
			p.push_back(exp1);
			cout << "\t\tExp1: " << exp1 << endl;
		}
		if (cfg.lookupValue("potential.exp2", exp2)) {
			p.push_back(exp2);
			cout << "\t\tExp2: " << exp2 << endl;
		}
		else {
			p.push_back(exp2);
			cout << "\t\tExp2: " << exp2 << endl;
		}
	}

	cout << endl;

	if (cfg.lookupValue("opt.name", algo)) {
		if (algo == "BFGS") {
			algo_switch = 1;
		}
		if (algo == "CG") {
			algo_switch = 2;
		}
		cout << "\tAlgo: " << algo << endl;
	}
	else {
	    cout << "\tNo 'opt' setting in configuration file." << endl;
		return 1;
	}

	if (cfg.lookupValue("opt.accuracy", accuracy)) {
		cout << "\t\tAccuracy: " << accuracy << endl;
	}
	else {
		cout << "\t\tAccuracy: " << accuracy << endl;
	}

	if (cfg.lookupValue("opt.dforce", dforce)) {
		cout << "\t\tDforce:  " << dforce << endl;
	}
	else {
		cout << "\t\tDforce:  " << dforce << endl;
	}

	if (cfg.lookupValue("opt.stepsize", stepsize)) {
		cout << "\t\tStepsize: " << stepsize << endl;
	}
	else {
		cout << "\t\tStepsize: " << stepsize << endl;
	}

	if (cfg.lookupValue("opt.nsteps", nsteps)) {
		cout << "\t\tNsteps: " << nsteps << endl;
	}
	else {
		cout << "\t\tNsteps: " << nsteps << endl;
	}
    opt.push_back(accuracy); //opt[0] should not be touched, not very important for opt
	opt.push_back(dforce);   //opt[1] for some reason can't be set below 10e-5, don't know why, GSL error
	opt.push_back(stepsize); //opt[2]
	opt.push_back(nsteps);   //opt[3]

	switches.push_back(potential_switch); 
	switches.push_back(algo_switch);
	switches.push_back(scaling_switch);
	switches.push_back(output_switch);

	cout << endl;

	return 0;
}
