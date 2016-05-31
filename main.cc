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

		//Inertia
		coord3d CoM = optKS[i].centreOfMass();
		optKS[i].shiftToCoM(CoM);
		vector< vector<double> > inertiaTensor = optKS[i].momentOfInertia();
		vector<double> inertia = diag(inertiaTensor);
		optKS[i].setMomentOfInertia(inertia);

		threadstream << "Eigenvalues of the hessian are:" << endl << eigenValues << endl;
		if (!optKS[i].isMinimum()) {
			threadstream << "Warning!!! Eigenvalue smaller than 0 in Hessian." << endl;
			hessianWarnings += 1;
			notMinimum.push_back(optKS[i].getNumber());
		}
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

	//SORT BY ENERGY
	sort(optKS.begin(), optKS.end());
	ofstream energies, energystats;

	energies.open ("energies");
	energystats.open("energystats");
	energies << left << setw(10) << "number" << right << setw(15) << "energy" << right << setw(25) << "eigenvalues" << endl;
    for (vector<structure>::size_type i = 0; i < optKS.size(); i++) {
		energies << 
			left << setw(10) << optKS[i].getNumber() << 
			right << setw(15) << optKS[i].getEnergy() << 
			right << setw(15) << optKS[i].getMomentOfInertia();
			if (find (notMinimum.begin(), notMinimum.end(), optKS[i].getNumber()) != notMinimum.end()) energies << "!!!";
		energies << endl;
	}

	//STATISTICS ON ENERGIES
	auto compare_inertia = [&] (vector<double> a, vector<double> b) {
		double threshold = 1e-5;
		if (b[0]-a[0] > threshold) return true;
		if (b[1]-a[1] > threshold) return true;
		if (b[2]-a[2] > threshold) return true;
		return false;
	};
	auto compare_map = [&] (pair < double, vector<double> > a, pair < double, vector<double> > b) { 
		return (b.first-a.first > 1e-10) && compare_inertia(a.second, b.second);
	};
	energyMap energyStat(compare_map);


	for (vector<structure>::size_type i = 0; i < optKS.size(); i++) {
		pair < double, vector<double> > key (optKS[i].getEnergy(), optKS[i].getMomentOfInertia());
		energyMap::iterator iter(energyStat.find(key));
		if (iter != energyStat.end()) {
			iter->second++;
		} else {
			energyStat[key] = 1;
		}	
	}


	energystats << "Statistics on energies: " << endl;
	energystats << setw(10) << right << "energy" << setw(30) << "count" << setw(10) << "inertia" << endl;
	for (energyMap::iterator iter = energyStat.begin(); iter != energyStat.end(); iter++) {
		energystats << setw(10) << right << iter->first.first << setw(30) << iter->second << setw(10) << iter->first.second << endl;
	}

	energystats << endl << "Number of unique energies is: " << energyStat.size() << endl;
	energystats << endl << "non-minimum structures:" << endl;
	for (vector<int>::size_type i = 0; i < notMinimum.size(); i++) {
		vector<structure>::iterator printThis = find_if (optKS.begin(), optKS.end(), [&] (structure toPrint) { return (toPrint.getNumber() == notMinimum[i]); });
		energystats << setw(10) << notMinimum[i] << setw(10) << printThis->getEnergy() << setw(10) << printThis->getMomentOfInertia();
		for (int j=0; j<10; j++) {
			energystats << ", " << printThis->getHessian()[j];
		}
		energystats << "..." << endl;
	}	

	energies.close();
	energystats.close();


	////WRITE OUTPUT STRUCTURES
	//switch (output_switch) {
	//	case 1:
	//		if (structureNumbers[0] == 0) {
	//			for (vector<structure>::size_type i = 0; i < optKS.size(); i++) {
	//				xyzout (optKS[i], "optStructure" + to_string (optKS[i].getNumber()));
	//				xyzout (allKS[i], "inpStructure" + to_string (allKS[i].getNumber()));
	//			}
	//		}
	//		else {
	//			for (vector<int>::size_type i = 0; i < structureNumbers.size(); i++) {
	//				vector<structure>::iterator printThis = find_if (optKS.begin(), optKS.end(), [&] (structure toPrint) { return (toPrint.getNumber() == structureNumbers[i]); });
	//				if (printThis != optKS.end()) {
	//					xyzout (*printThis, "optStructure" + to_string (printThis->getNumber()));
	//					xyzout (allKS[structureNumbers[i - 1]], "inpStructure" + to_string (allKS[i].getNumber()));
	//				}
	//				else {
	//					cerr << "Structure number " << structureNumbers[i] << " not found." << endl;
	//				}
	//			}
	//		}
	//	default:
	//		break;
	//}

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
