#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <functional>
#include <cassert>
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






	//CHECK IF FILE EXISTS
    string fileName = "coord"; //safe input file name
    
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
    vector<structure> optKS = readallstruct(fileName);
	
	//if scaling is found in settings file scale all coordinates accordingly
	switch (scaling_switch) {
		case 1:
			for (vector<structure>::size_type i = 0; i < optKS.size(); i++) {
				optKS[i] *= scalingFactor;
			}
		default:
			break;
	}






	//HESSIAN
	ofstream energystats;
	unsigned int hessianWarnings = 0;
	vector<int> notMinimum;
	vector<structure> notMinimumKS;

	#pragma omp parallel for
	for (vector<structure>::size_type i = 0; i < optKS.size(); i++) {

		vector< vector<double> > hessian;
		vector<double> eigenValues;
		stringstream threadstream;

		optKS[i].setEnergy (optKS[i].sumOverAllInteractions(p));

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
		#pragma omp critical
		{
			if (!optKS[i].isMinimum()) {
				threadstream << "Warning!!! Eigenvalue smaller than 0 in Hessian." << endl;
				hessianWarnings += 1;
				notMinimum.push_back(optKS[i].getNumber());
				notMinimumKS.push_back(optKS[i]);
			}
		energystats << threadstream.rdbuf() << endl;
		}
	}




	topt=clock();
	float optTime ((float)topt-(float)tstart);
	cout << "\tTime for Hessian and energy calculation: " << optTime/CLOCKS_PER_SEC << " s" << endl << endl;


	





	//SORT BY ENERGY
	sort(optKS.begin(), optKS.end());
	ofstream energies;

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
	auto compare_map = [&] (pair < double, vector<double> > a, pair < double, vector<double> > b) { //sort by energy and inertia function
		double eps = 1e-3;
		if (b.first-a.first > eps) return true;
		if (a.first-b.first < eps && b.second[0]-a.second[0] > eps) return true;
		if (a.first-b.first < eps && a.second[0]-b.second[0] < eps && b.second[1]-a.second[1] > eps) return true;
		if (a.first-b.first < eps && a.second[0]-b.second[0] < eps && a.second[1]-b.second[1] < eps && b.second[2]-a.second[2] > eps) return true;	
		return false;
	};

	energyMap energyStat(compare_map);
	energyMap notMinimumStat(compare_map);

	//put all structures in a map with energy and inertia as key
	for (vector<structure>::size_type i = 0; i < optKS.size(); i++) {
		pair < double, vector<double> > key (optKS[i].getEnergy(), optKS[i].getMomentOfInertia());
		energyMap::iterator iter(energyStat.find(key));
		if (iter != energyStat.end()) {
			iter->second++;
		} else {
			energyStat[key] = 1;
		}	
	}

	//put all nonMinimum structures in a map with energy and inertia as key
	for (vector<structure>::size_type i = 0; i < notMinimumKS.size(); i++) {
		pair < double, vector<double> > key (notMinimumKS[i].getEnergy(), notMinimumKS[i].getMomentOfInertia());
		energyMap::iterator iter(notMinimumStat.find(key));
		if (iter != notMinimumStat.end()) {
			iter->second++;
		} else {
			notMinimumStat[key] = 1;
		}	
	}


	energystats << "Statistics on energies: " << endl;
	energystats << setw(10) << right << "energy" << setw(30) << "count" << setw(10) << "inertia" << endl;
	for (energyMap::iterator iter = energyStat.begin(); iter != energyStat.end(); iter++) {
		energystats << setw(10) << right << iter->first.first << setw(30) << iter->second << setw(10) << iter->first.second << endl;
	}

	energystats << endl << "non-minimum structures:" << endl;
	for (energyMap::iterator iter = notMinimumStat.begin(); iter != notMinimumStat.end(); iter++) {
		energystats << setw(10) << right << iter->first.first << setw(30) << iter->second << setw(10) << iter->first.second << endl;
	}
	energystats << endl << "Number of unique structures: " << energyStat.size() << endl;
	energystats << "Number of unique non-minimum structures: " << notMinimumStat.size() << endl;
	energystats << "Number of unique minimum structures: " << energyStat.size()-notMinimumStat.size() << endl;
	//for (vector<int>::size_type i = 0; i < notMinimum.size(); i++) {
	//	vector<structure>::iterator printThis = find_if (optKS.begin(), optKS.end(), [&] (structure toPrint) { return (toPrint.getNumber() == notMinimum[i]); });
	//	energystats << setw(10) << notMinimum[i] << setw(10) << printThis->getEnergy() << setw(10) << printThis->getMomentOfInertia();
	//	for (int j=0; j<10; j++) {
	//		energystats << ", " << printThis->getHessian()[j];
	//	}
	//	energystats << "..." << endl;
	//}	

	energystats << endl << "Number of warnings for Hessian eigenvalues: " << hessianWarnings << endl;
	energystats << "List of non-minimum Structures: " << notMinimum << endl;

	energies.close();
	energystats.close();

	//test for inversion matrix
	matrix3d test = optKS[0].m3d_momentOfInertia();
	matrix3d test1 = m3d_diagv(test);
	matrix3d test2 = test1.inverse();

	for (int i=0; i<3; i++) 
		for (int j=0; j<3; j++)
			cout << i << j << " " << test2(i,j) << endl;


	//calculate all interparticle distances for minimum structures
	for (vector<structure>::size_type i=0; i<optKS.size(); i++) {
		vector<double> interPartDist;
		vector<coord3d> currentCoord=optKS[i].getCoordinates();
		for (vector<coord3d>::size_type j=0; j<currentCoord.size(); j++) {
			for (vector<coord3d>::size_type k=j+1; k<currentCoord.size(); k++) {
				interPartDist.push_back(coord3d::dist(currentCoord[j], currentCoord[k]));	
			}
		}
		sort(interPartDist.begin(), interPartDist.end());
		optKS[i].setInterPartDist(interPartDist);
	}
	
	//calculate all interparticle distances for non-minimum structures
	for (vector<structure>::size_type i=0; i<notMinimumKS.size(); i++) {
		vector<double> interPartDist;
		vector<coord3d> currentCoord=notMinimumKS[i].getCoordinates();
		for (vector<coord3d>::size_type j=0; j<currentCoord.size(); j++) {
			for (vector<coord3d>::size_type k=j+1; k<currentCoord.size(); k++) {
				interPartDist.push_back(coord3d::dist(currentCoord[j], currentCoord[k]));	
			}
		}
		sort(interPartDist.begin(), interPartDist.end());
		notMinimumKS[i].setInterPartDist(interPartDist);
	}

	//compare function for distance vectors
	auto compare_interPartDist = [&] (vector<double> a, vector<double> b) {
		assert (a.size() == b.size());
		double eps = 1e-3;
		vector<double> diff;
		for (vector<double>::size_type i = 0; i < a.size(); i++) {
			diff.push_back(abs(a[i]-b[i]));
		}
		auto max = max_element (begin(diff), end(diff));
		if (*max < eps) return true;
		return false;
	};

	//sort all minimum structures into vector< vector <structure> > according to distance vector
	vector< vector<structure> > eqClasses;
	vector<structure> one;
	one.push_back(optKS[0]);
	eqClasses.push_back(one);
	for (vector<structure>::size_type i = 1; i < optKS.size(); i++) {
		for (vector< vector<structure> >::size_type j = 0; j < eqClasses.size(); j++) {
			if (compare_interPartDist (optKS[i].getInterPartDist(), eqClasses[j][0].getInterPartDist())) {
				eqClasses[j].push_back(optKS[i]);
				break;
			}
			if (j + 1 == eqClasses.size()) {
				vector<structure> newEqClass;
				newEqClass.push_back(optKS[i]);
				eqClasses.push_back(newEqClass);
			}
		}
	}
	cout << "Number of equality classes: " << eqClasses.size() << endl;


	//sort all non-minimum structures into vector< vector <structure> > according to distance vector
	vector< vector<structure> > notMinimumEqClasses;
	vector<structure> notMinimumOne;
	notMinimumOne.push_back(notMinimumKS[0]);
	notMinimumEqClasses.push_back(notMinimumOne);
	for (vector<structure>::size_type i = 1; i < notMinimumKS.size(); i++) {
		for (vector< vector<structure> >::size_type j = 0; j < notMinimumEqClasses.size(); j++) {
			if (compare_interPartDist (notMinimumKS[i].getInterPartDist(), notMinimumEqClasses[j][0].getInterPartDist())) {
				notMinimumEqClasses[j].push_back(notMinimumKS[i]);
				break;
			}
			if (j + 1 == notMinimumEqClasses.size()) {
				vector<structure> newEqClass;
				newEqClass.push_back(notMinimumKS[i]);
				notMinimumEqClasses.push_back(newEqClass);
			}
		}
	}
	cout << "Number of non-minimum equality classes: " << notMinimumEqClasses.size() << endl;

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





	
	
	tend=clock();
	float totalTime ((float)tend-(float)tstart);
	cout << "\tTotal rumtime: " << totalTime/CLOCKS_PER_SEC << " s" << endl;
    return 0; 

    
}
