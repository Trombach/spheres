#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include "spheres.h"
#include <libconfig.h++>

using namespace std; 
using namespace libconfig;
//
//function for determining if a line in the input is emtpy
//used for breaking the read loop
//
bool justempty(string str) {
	if (str.empty()) {
		return true;
	}
	return false;
}

//
//function for testing if input file exists
//returns true for existing file
//
bool fexists (const std::string& fileName) {
    ifstream infile(fileName.c_str());
    return infile;
}

//
//function for reading all strucutures of input file
//structures must be separated by a blank line, last line should be a blank line as well
//
vector<structure> readallstruct (const std::string& fileName) {
	ifstream infile(fileName.c_str());
	coord3d sphere; //sphere refers to 1 'atom'
	string line;
	vector<structure> allKissingSpheres; //this vector stores all read in structures
	//double while loop, inner loop goes over all coordinates of one structure, outer loop goes over all 
	//structures till eof
	while (true) {
	    structure kissingSphere; //kissingSphere refers to 1 cluster of kissing Spheres
	    while (true) {
	        getline(infile,line);
	        //cout << line << endl;
	        if(justempty(line) || infile.eof()) break;
	        else {
	    	    stringstream lineStream(line); //string can't be read into coord3d object, convert to 
				                               //stringstream first
	    	    lineStream >> sphere;
	            kissingSphere.push_back(sphere);
	            //cout << sphere << endl;
	        }
	    }
		if (infile.eof()) break; //does last line has to be blank? break works for files with blank line at 
		                         //the end
		else {
		    allKissingSpheres.push_back(kissingSphere);
            //cout << "Strucuture has been read." << endl;
	    }
	}
    infile.close();
	cout << "Number of structures is " << allKissingSpheres.size() << endl;
	//cout << allKissingSpheres[0][3] << endl;
	//cout << allKissingSpheres[1][3] << endl;
	cout << "Last structure is:" << endl;
	for (structure::size_type i = 0; i<allKissingSpheres [allKissingSpheres.size() - 1].size(); ++i) {
	    cout <<  allKissingSpheres [allKissingSpheres.size() - 1][i]  << endl;
	}
	return allKissingSpheres;
}

//
//function for reading and processing the settings file
//
int readSettings () {
    ifstream infile("settings");
	string line;
	
	while (infile >> line) {
        cout << line << endl;
		string potential;
	    if (line == "%pot") {
		    getline (infile, line);
			potential = line;
		    cout << "potential is: " << potential << endl;

			if (potential == "LJ") {
				double epsilon;
				//double rm;
			    while (infile >> line) {
                    if (line == "epsilon") {
				        getline (infile, line);
						epsilon = std::stod(line);
						cout << "Epsilon is: " << epsilon << endl; 
			        }
				}	
		    }	
			else cout << "Potential unknown" << endl;
			return 1;
		}	
	    	
    }
	return 0;
}	
//
//MAIN FUNCTION BEGINS HERE
//

int main (int argc, char *argv[]) {
	
//
//SOME START UP CHECKS AND ARGUMENT PROCESSING
//
    if ( argc < 2 ) {
        cout << "Please provide filename." << endl;
        return 1;
    }
    else {
        cout << "Filename is " << argv[ argc -1 ] << "."  << endl;
        //cout << "Number of arguments is " << argc << endl;
    }

//
//CHECK IF FILE EXISTS
//

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

//
//READ IN ALL STRUCTURES AT ONCE AND CALCULATE LJ-ENERGY FOR EACH STRUCTURE
//
    vector<structure> allKissingSpheres = readallstruct(fileName);
	vector<double> allEnergies;
    for (vector<structure>::size_type i = 0; i < allKissingSpheres.size(); ++i) {
	    allEnergies.push_back( allKissingSpheres[i].sumOverAllInteractions() );
	}

    cout << "Number of Energies is " << allEnergies.size() << endl;	
	for (vector<double>::size_type i=0; i < 1; ++i) {
	    cout << "Total LJ-Energy for structure " << i + 1 << " is " << allEnergies[i] << endl;
    }
    
	cout << "Structure is:" << endl;
	for (structure::size_type i = 0; i < allKissingSpheres[0].size(); ++i) {
		cout << allKissingSpheres[0][i] << endl;
	}
//
//CALCULATE GRADIENTS
//
//	vector< vector<coord3d> > allGradients;
//    for (vector<structure>::size_type i = 0; i < allKissingSpheres.size(); ++i) {
//        vector<coord3d> gradients = allKissingSpheres[i].sumOverAllGradients();
//        allGradients.push_back(gradients);
//	}
//    for (vector< vector<coord3d> >::size_type i = 0; i < allGradients.size(); ++i) { 
//        coord3d sum(0,0,0);
//	    for (structure::size_type j = 0; j < allGradients[i].size(); ++j) {
//	    	sum += allGradients[i][j];
//	    }
//	    cout << "Sum over all Forces for Structure " << i + 1 << " is : " << sum << endl;
//	}
//
//
//READ SETTINGS FILE
//

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
	string potential;
	vector<double> p;
	int potential_switch;
	int algo_switch;
    string algo;
	double accuracy;
	double dforce;
	double stepsize;
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
		double epsilon;
		double rm;
		cout << "Setting LJ Parameters." << endl;
		if (cfg.lookupValue("potential.epsilon", epsilon)) {
			p.push_back(epsilon);
			cout << "Epsilon set to " << p[0] << endl;
		}
		else {
			cout << "No 'epsilon' in setting file." << endl;
			return 1;
		}
		if (cfg.lookupValue("potential.rm", rm)) {
			p.push_back(rm);
			cout << "Rm set to " << p[1] << endl;
		}
		else {
			cout << "No 'rm' in setting file." << endl;
			return 1;
		}
	}
	vector<double> opt; //vector of algo settings, 0 == accuracy, 1 == dforce, 2 == stepsize
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
		accuracy = 1e-4;
		cout << "Fall back to " << accuracy << endl;
	}

	if (cfg.lookupValue("opt.dforce", dforce)) {
		cout << "Convergence criterion for gradients set to " << dforce << endl;
	}
	else {
	    cout << "No 'dforce' setting in configuration file." << endl;
		dforce = 1e-3;
	    cout << "Fall back to " << dforce << endl;
	}

	if (cfg.lookupValue("opt.stepsize", stepsize)) {
		cout << "Stepsize set to " << stepsize << endl;
	}
	else {
	    cout << "No 'stepsize' setting in configuration file." << endl;
		stepsize = 0.01;
		cout << "Fall back to " << stepsize << endl;
	}
    opt.push_back(accuracy);
	opt.push_back(dforce);
	opt.push_back(stepsize);
    

    allKissingSpheres[0].optimize(algo_switch, potential_switch, p, opt);    

    return 0; 

    
}