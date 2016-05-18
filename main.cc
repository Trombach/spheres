#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <libconfig.h++>
#include <algorithm>
#include <map>
#include "structure.h"
#include "iop.h"
#include "lina.h"


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

typedef map<double, unsigned int, function<bool(double a, double b)> > energyMap;

//MAIN FUNCTION BEGINS HERE

int main (int argc, char *argv[]) {
	
	cout << endl;
	//SOME START UP CHECKS AND ARGUMENT PROCESSING
    if ( argc < 2 ) {
        cout << "\t" << "Please provide filename!" << endl;
        return 1;
    }

	//CHECK IF FILE EXISTS

    string fileName = argv[ argc -1 ]; //safe input file name
    
    if (fexists(fileName)) {
        cout << "\t" << "File " << fileName << " exists." << endl;
    }
    else {
        cout << "\t" << "File " << fileName << " does not exist." << endl;
        return 1;
    }

	cout << endl;

    if (fexists("settings")) {
		cout << "\t" << "Settings file found." << endl;
	}
	else {
		cout << "\t" << "No settings file in working directory" << endl;
		return 1;
	}

	
	//READ SETTINGS FILE

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

	cout << endl;

	const Setting &root = cfg.getRoot();

	string potential;
	vector<double> p;
	vector<int> structureNumbers;
	int potential_switch, algo_switch, scaling_switch, output_switch(1);
    string algo;
	double scalingFactor(1.0), accuracy(0.1), dforce(1e-3), stepsize(0.01);
	int nsteps(100);
	bool printInput(false);
	
	try {
		const Setting &numbers = root["output"]["number"];
		int count = numbers.getLength();
		for (int i = 0; i < count; i++) {
			structureNumbers.push_back (numbers[i]);
		}
	}
	catch (const SettingNotFoundException &nfex) {
		cerr << "\tOutput setting not Found." << endl;
		output_switch = 0;
	}
	cout << "\tStructures to print: " << structureNumbers;
	if (cfg.lookupValue ("output.input", printInput)) {
		cout << ", including input coordinates" << endl;
	}
	else {
		cout << ", excluding input coordinates" << endl;
	}


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

	vector<double> opt; //vector of algo settings, 0 == accuracy, 1 == dforce, 2 == stepsize, 3 == nsteps
	if (cfg.lookupValue("opt.name", algo)) {
		if (algo == "BFGS") {
			algo_switch = 1;
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

	cout << endl;
    
	//READ IN ALL STRUCTURES AT ONCE
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


	//OPTIMIZE AND HESSIAN
	cout << "\tStarting optimization..." << endl;
	ofstream min;
	unsigned int hessianWarnings = 0;
	vector<int> notMinimum;
	min.open ("opt");
	for (vector<structure>::size_type i = 0; i < allKissingSpheres.size(); i++) {
		min << "Optimization for structure no " << allKissingSpheres[i].getNumber() << endl;
        optimizedKissingSpheres.push_back (allKissingSpheres[i].optimize(min, algo_switch, potential_switch, p, opt , allEnergies));
		optimizedKissingSpheres[i].setNumber( allKissingSpheres[i].getNumber() );
		hessian = optimizedKissingSpheres[i].hessian(p);
		eigenValues = diag(hessian);
		optimizedKissingSpheres[i].setHessian (eigenValues);
		min << "Eigenvalues of the hessian are:" << endl << eigenValues << endl;
		if (!optimizedKissingSpheres[i].isMinimum()) {
			min << "Warning!!! Eigenvalue smaller than 0 in Hessian." << endl;
			hessianWarnings += 1;
			notMinimum.push_back(optimizedKissingSpheres[i].getNumber());
		}
		min << "###############################################################\n" << endl;
	}
	min << endl << "Number of warnings for Hessian eigenvalues: " << hessianWarnings << endl;
	min << "List of non-minimum Structures: " << notMinimum << endl;
	min.close();
	cout << "\tAll Done!" << endl << endl;
	
	//INERTIA TENSOR AND EIGENVALUES
	for (vector<structure>::size_type i = 0; i < optimizedKissingSpheres.size(); i++) {
		coord3d CoM = optimizedKissingSpheres[i].centreOfMass();
		optimizedKissingSpheres[i].shiftToCoM(CoM);
		vector< vector<double> > inertiaTensor = optimizedKissingSpheres[i].momentOfInertia();
		vector<double> inertia = diag(inertiaTensor);
		optimizedKissingSpheres[i].setMomentOfInertia(inertia);
	}

	//SORT BY ENERGY
	sort(optimizedKissingSpheres.begin(), optimizedKissingSpheres.end());
	ofstream energies, energystats;

	auto compare_map = [&] (double a, double b) { return b-a > 0.00000001;};
	energyMap energyStat(compare_map);

	energies.open ("energies");
	energystats.open("energystats");
	energies << left << setw(10) << "number" << right << setw(15) << "energy" << right << setw(25) << "eigenvalues" << endl;
    for (vector<structure>::size_type i = 0; i < optimizedKissingSpheres.size(); i++) {
		energies << 
			left << setw(10) << optimizedKissingSpheres[i].getNumber() << 
			right << setw(15) << optimizedKissingSpheres[i].getEnergy() << 
			right << setw(15) << optimizedKissingSpheres[i].getMomentOfInertia();
			if (find (notMinimum.begin(), notMinimum.end(), optimizedKissingSpheres[i].getNumber()) != notMinimum.end()) energies << "!!!";
		energies << endl;
	}


	//STATISTICS ON ENERGIES
	for (vector<structure>::size_type i = 0; i < optimizedKissingSpheres.size(); i++) {
		double energy = optimizedKissingSpheres[i].getEnergy();
		energyMap::iterator iter(energyStat.find(energy));
		if (iter != energyStat.end()) {
			iter->second++;
		} else {
			energyStat[energy] = 1;
		}	
	}
	energystats << "Statistics on energies: " << endl;
	energystats << setw(10) << "energy" << setw(10) << "count" << endl;
	for (energyMap::iterator iter = energyStat.begin(); iter != energyStat.end(); iter++) {
		energystats << setw(10) << iter->first << setw(10) << iter->second << endl;
	}
	energystats << endl << "Number of unique energies is: " << energyStat.size() << endl;
	energystats << endl << "non-minimum structures:" << endl;
	for (vector<int>::size_type i = 0; i < notMinimum.size(); i++) {
		vector<structure>::iterator printThis = find_if (optimizedKissingSpheres.begin(), optimizedKissingSpheres.end(), [&] (structure toPrint) { return (toPrint.getNumber() == notMinimum[i]); });
		energystats << setw(10) << notMinimum[i] << setw(10) << printThis->getEnergy() << setw(10) << printThis->getMomentOfInertia() << endl;
	}	
	energies.close();
	energystats.close();


	//WRITE OUTPUT STRUCTURES
	switch (output_switch) {
		case 1:
			if (structureNumbers[0] == 0) {
				for (vector<structure>::size_type i = 0; i < optimizedKissingSpheres.size(); i++) {
					xyzout (optimizedKissingSpheres[i], "optStructure" + to_string (optimizedKissingSpheres[i].getNumber()));
					xyzout (allKissingSpheres[i], "inpStructure" + to_string (allKissingSpheres[i].getNumber()));
				}
			}
			else {
				for (vector<int>::size_type i = 0; i < structureNumbers.size(); i++) {
					vector<structure>::iterator printThis = find_if (optimizedKissingSpheres.begin(), optimizedKissingSpheres.end(), [&] (structure toPrint) { return (toPrint.getNumber() == structureNumbers[i]); });
					if (printThis != optimizedKissingSpheres.end()) {
						xyzout (*printThis, "optStructure" + to_string (printThis->getNumber()));
						xyzout (allKissingSpheres[structureNumbers[i - 1]], "inpStructure" + to_string (allKissingSpheres[i].getNumber()));
					}
					else {
						cerr << "Structure number " << structureNumbers[i] << " not found." << endl;
					}
				}
			}
		default:
			break;
	}

	cout << "\tAll good!" << endl;
    return 0; 

    
}
