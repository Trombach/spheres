#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include "spheres.h"


using namespace std; 

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
    
//
//READ IN ALL STRUCTURES AT ONCE AND CALCULATE LJ-ENERGY FOR EACH STRUCTURE
//
    vector<structure> allKissingSpheres = readallstruct(fileName);
	vector<double> allEnergies;
    for (vector<structure>::size_type i = 0; i < allKissingSpheres.size(); ++i) {
	    allEnergies.push_back( allKissingSpheres[i].sumOverAllInteractions() );
	}

    cout << "Number of Energies is " << allEnergies.size() << endl;	
	for (vector<double>::size_type i=0; i < allEnergies.size(); ++i) {
	    cout << "Total LJ-Energy for structure " << i + 1 << " is " << allEnergies[i] << endl;
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
//OPTIMIZE STRUCTURE
//
    allKissingSpheres[0].optimize();    

    return 0;  

}
