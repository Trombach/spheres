#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include "spheres.h"


using namespace std;
string fileName;

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
	//double while loop, inner loop goes over all coordinates of one structure, outer loop goes over all structures till eof
	while (true) {
	    structure kissingSphere; //kissingSphere refers to 1 cluster of kissing Spheres
	    while (true) {
	        getline(infile,line);
	        //cout << line << endl;
	        if(justempty(line) || infile.eof()) break;
	        else {
	    	    stringstream lineStream(line); //string can't be read into coord3d object, convert to stringstream first
	    	    lineStream >> sphere;
	            kissingSphere.push_back(sphere);
	            //cout << sphere << endl;
	        }
	    }
		if (infile.eof()) break; //does last line has to be blank? break works for files with blank line at the end
		else {
		    allKissingSpheres.push_back(kissingSphere);
            //cout << "Strucuture has been read." << endl;
	    }
	}

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
//function to calculate LJ Energy, needs rij,rm and epsilon
//
double LJEnergy (const double distance,const double epsilon,const double rm) {
    return epsilon * ( (pow (rm / distance, 12)) - 2 * (pow (rm / distance, 6)) );
}

//
//function to sum over all sphere interactions, change later to work with different potentials
//
double sumOverAllInteractions (structure kissingSphere) {
    double totalEnergy = 0;
	for (structure::iterator iter = kissingSphere.begin(); iter != kissingSphere.end(); ++iter) { 
		for (structure::iterator jter = iter + 1; jter != kissingSphere.end(); ++jter) {
			//cout << "iter is " << *iter << endl;
			//cout << "jter is " << *jter << endl;
            totalEnergy += LJEnergy (coord3d::dist (*iter,*jter),1,0.5);
        }
	}
	return totalEnergy;
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

    fileName = argv[ argc -1 ]; //safe input file name
    cout << fileName << endl;
    
    if (fexists(fileName)) {
        cout << "File " << fileName << " exists." << endl;
    }
    else {
        cout << "File " << fileName << " does not exist." << endl;
        return 1;
    }
    
//
//READ IN ALL STRUCTURES AT ONCE AND CALCULATE LJ-ENERGY FOR EACH STRUCTURE
//
    //iterate over double index ij, where N>j>i and N>i
    vector<structure> allKissingSpheres = readallstruct(fileName);
	vector<double> allEnergies;
    for (vector<structure>::size_type i = 0; i < allKissingSpheres.size(); ++i) {
	    allEnergies.push_back(sumOverAllInteractions(allKissingSpheres[i]));
	}

    cout << "Number of Energies is " << allEnergies.size() << endl;	
	for (vector<double>::size_type i=0; i < allEnergies.size(); ++i) {
	    cout << "Total LJ-Energy or structure " << i + 1 << " is " << allEnergies[i] << endl;
    }
    return 0;  

}
