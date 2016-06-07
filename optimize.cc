#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include "structure.h"
#include "iop.h"
#include "lina.h"


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

typedef map <pair < double, vector<double> >, unsigned int, function<bool( pair < double, vector<double> > a, pair < double, vector<double> > b)> > energyMap;

//MAIN FUNCTION BEGINS HERE

int main (int argc, char *argv[]) {
	clock_t tstart, tend, topt;
	tstart=clock();




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






	//READ SETTINGS FILE
    if (fexists("settings")) {
		cout << "\t" << "Settings file found." << endl;
	}
	else {
		cout << "\t" << "No settings file in working directory" << endl;
		return 1;
	}

	vector<double> p;
	vector<double> opt; //vector of algo settings, 0 == accuracy, 1 == dforce, 2 == stepsize, 3 == nsteps
	vector<int> switches; //vector of switches, 0 == potential, 1 == algo, 2 == scaling
	double scalingFactor(1.0);

	readsettings(opt, p, switches, scalingFactor);
	int potential_switch(switches[0]), algo_switch(switches[1]), scaling_switch(switches[2]);
	
	
	cout << endl;






    
	//READ IN ALL STRUCTURES AT ONCE
    vector<structure> allKS = readallstruct(fileName);
	
	//if scaling is found in settings file scale all coordinates accordingly
	switch (scaling_switch) {
		case 1:
			for (vector<structure>::size_type i = 0; i < allKS.size(); i++) {
				allKS[i] *= scalingFactor;
			}
		default:
			break;
	}






	//OPTIMIZE AND HESSIAN
    vector<structure> optKS; optKS.resize(allKS.size());


	cout << "\tStarting optimization..." << endl;
	ofstream min;
	unsigned int hessianWarnings = 0;
	vector<int> notMinimum;
	min.open ("opt");

	#pragma omp parallel for
	for (vector<structure>::size_type i = 0; i < allKS.size(); i++) {

		vector< vector<double> > hessian;
		vector<double> eigenValues;
		stringstream threadstream;

		threadstream << "Optimization for structure no " << allKS[i].getNumber() << endl;

        optKS[i] = allKS[i].optimize(threadstream, algo_switch, potential_switch, p, opt);
		optKS[i].setNumber( allKS[i].getNumber() );

		hessian = optKS[i].hessian(p);
		eigenValues = diag(hessian);
		optKS[i].setHessian (eigenValues);


		threadstream << "Eigenvalues of the hessian are:" << endl << eigenValues << endl;

		int reopts(0);
		vector<coord3d> gradients;
		vector<coord3d> coords;
		while (!optKS[i].isMinimum() && reopts < 6) {
			reopts++;
			threadstream << "Reoptimization attempt " << reopts << endl;
			gradients = optKS[i].sumOverAllGradients(p);
			coords = optKS[i].getCoordinates();

			if (gradients.size() != coords.size()) {
				cerr << "gradient and coordinate vector are of different size" << endl;
				break;
			}

			for (vector<coord3d>::size_type	j = 0; j < coords.size(); j++) {
				coords[j] = coords[j] + gradients[j];
			}
			optKS[i].setCoordinates(coords);
			
			optKS[i] = optKS[i].optimize(threadstream, algo_switch, potential_switch, p, opt);

			hessian = optKS[i].hessian(p);
			eigenValues = diag(hessian);
			optKS[i].setHessian (eigenValues);
		}
		if (!optKS[i].isMinimum() || reopts > 5) {
			threadstream << "Structure did not converge to minimum after " << reopts << " attempts." << endl;
			threadstream << "Final gradient:" << endl;
			for (vector<coord3d>::size_type s = 0; s < gradients.size(); s++) {
				threadstream << gradients[s] << endl;
			}
			hessianWarnings +=1;
			notMinimum.push_back(optKS[i].getNumber());
		}

		//Inertia
		coord3d CoM = optKS[i].centreOfMass();
		optKS[i].shiftToCoM(CoM);
		vector< vector<double> > inertiaTensor = optKS[i].momentOfInertia();
		vector<double> inertia = diag(inertiaTensor);
		optKS[i].setMomentOfInertia(inertia);

		#pragma omp critical
		{
			min << threadstream.rdbuf() << endl;
			min << "***************End of Opt***************" << endl;
		}
	}
	





	min << endl << "Number of warnings for Hessian eigenvalues: " << hessianWarnings << endl;
	min << "List of non-minimum Structures: " << notMinimum << endl;
	min.close();






	topt=clock();
	float optTime ((float)topt-(float)tstart);
	cout << "\tTime for structure optimization: " << optTime/CLOCKS_PER_SEC << " s" << endl << endl;






	//
	//WRITE COORD CHECKPOINT
	
	stringstream out;
	simpleout(optKS, out);

	ofstream coord;
	coord.open("coord");
	coord << out.rdbuf();
	coord.close();
	
	
	tend=clock();
	float totalTime ((float)tend-(float)tstart);
	cout << "\tTotal rumtime: " << totalTime/CLOCKS_PER_SEC << " s" << endl;
    return 0; 

    
}
