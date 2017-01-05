#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <functional>
#include "structure.h"
#include "iop.h"
#include "lina.h"
#include "timer.h"


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


//MAIN FUNCTION BEGINS HERE

int main (int argc, char *argv[]) {

	timer T;


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
	cout << "\t#structures :" << allKS.size() << endl;
	
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
    vector<structure> optKS; //optKS.resize(allKS.size());


	cout << "\tStarting optimization..." << endl;
	ofstream min;
	unsigned int hessianWarnings = 0;
	vector<int> notMinimum;
	vector<structure> notMinimumKS;
	min.open ("opt");
	min.precision(14);

	#pragma omp parallel for
	for (vector<structure>::size_type i = 0; i < allKS.size(); i++) {

		vector< vector<double> > hessian;
		vector<double> eigenValues;
		structure threadKS;
		stringstream threadstream;

		threadstream.precision(14);

		threadstream << "Optimization for structure no " << allKS[i].getNumber() << endl;

        threadKS = allKS[i].optimize(threadstream, algo_switch, potential_switch, p, opt);
		threadKS.setNumber( allKS[i].getNumber() );


		hessian = threadKS.hessian(p);
		eigenValues = diag(hessian);
		threadKS.setHessian (eigenValues);


		//threadstream << "Eigenvalues of the hessian are:" << endl << eigenValues << endl;




		#pragma omp critical
		{
			if (threadKS.isMinimum()) {
				optKS.push_back(threadKS);
			}
			else {
				hessianWarnings += 1;
				notMinimum.push_back(threadKS.getNumber());
				notMinimumKS.push_back(threadKS);
				threadstream << "Structure did not converge to minimum, will reoptimize later." << endl;
			}
			min << threadstream.rdbuf() << endl;
			min << "***************End of Opt***************" << endl;
		}
	}





	min << endl << "Number of warnings for Hessian eigenvalues: " << hessianWarnings << endl;
	min << "List of non-minimum Structures: " << notMinimum << endl;
	min.close();

	
	cout << "\tStarting reoptimization for " << notMinimumKS.size() << " structure(s)." << endl;
	//ATTEMPT REOPTIMIZATION
	//positive
	
	ofstream remin;
	remin.open("reopt");
	vector<structure>  reoptKS;

	remin << endl << "////POSITIVE////" << endl << endl;
	for (vector<structure>::size_type i=0; i< notMinimumKS.size(); i++) {
		structure KS = notMinimumKS[i], newKS;
		stringstream threadstream;
		int reopts(0);


		threadstream << "Reoptimization in positive direction for structure No. " << notMinimumKS[i].getNumber() << endl;


		do {
			reopts++;
			vector<coord3d> gradients;
			vector<coord3d> coords;
			vector<vector<double> > hessian;
			vector<double> eigenValues;
			

			threadstream << "Reoptimization attempt " << reopts << endl;
			coords = KS.getCoordinates();
			hessian = KS.hessian(p);
			vector<pair<double, vector<double> > > eval_evec = diagv(hessian);

			auto pairCompare = [&] (const pair<double, vector<double> > a, const pair<double, vector<double> > b) {
				return a.first < b.first;
			};
			sort (eval_evec.begin(), eval_evec.end(), pairCompare);
			
			for (vector<double>::size_type j=0; j<eval_evec[0].second.size()/3; j++) {
				coord3d gradient (eval_evec[0].second[3*j], eval_evec[0].second[3*j+1], eval_evec[0].second[3*j+2]);
				gradients.push_back(gradient / gradient.norm());
			}

			if (gradients.size() != coords.size()) {
				cerr << "gradient and coordinate vector are of different size" << endl;
				break;
			}

			vector<coord3d> displacement;

			threadstream << "negative deflection" << endl;
			for (vector<coord3d>::size_type	j = 0; j < coords.size(); j++) {
				displacement.push_back(coords[j] + (gradients[j] * 0.1));
			}

			KS.setCoordinates(displacement);

			newKS = KS.optimize(threadstream, algo_switch, potential_switch, p, opt);

			hessian = newKS.hessian(p);
			eigenValues = diag(hessian);
			newKS.setHessian (eigenValues);
			
			newKS.setNumber(notMinimumKS[i].getNumber());


			//Inertia
			coord3d CoM = newKS.centreOfMass();
			newKS.shiftToCoM(CoM);
			vector< vector<double> > inertiaTensor = newKS.momentOfInertia();
			vector<double> inertia = diag(inertiaTensor);
			newKS.setMomentOfInertia(inertia);

		}
		while (reopts < 6 && !newKS.isMinimum());
		
		if (!newKS.isMinimum()) {
			threadstream << "Structure did not converge after " << reopts << " attempts." << endl;
			cerr << "Reopt failure, structure " << newKS.getNumber() << ", positive displacement" << endl;
		}
		reoptKS.push_back(newKS);
		threadstream << "***************End of Opt***************" << endl;
		remin << threadstream.rdbuf() << endl;	
	}

	remin << endl << "////NEGATIVE////" << endl << endl;

	for (vector<structure>::size_type i=0; i< notMinimumKS.size(); i++) {
		structure KS = notMinimumKS[i], newKS;
		stringstream threadstream;
		int reopts(0);


		threadstream << "Reoptimization in negative direction for structure No. " << notMinimumKS[i].getNumber() << endl;


		do {
			reopts++;
			vector<coord3d> gradients;
			vector<coord3d> coords;
			vector<vector<double> > hessian;
			vector<double> eigenValues;
			

			threadstream << "Reoptimization attempt " << reopts << endl;
			coords = KS.getCoordinates();
			hessian = KS.hessian(p);
			vector<pair<double, vector<double> > > eval_evec = diagv(hessian);

			auto pairCompare = [&] (const pair<double, vector<double> > a, const pair<double, vector<double> > b) {
				return a.first < b.first;
			};
			sort (eval_evec.begin(), eval_evec.end(), pairCompare);
			
			for (vector<double>::size_type j=0; j<eval_evec[0].second.size()/3; j++) {
				coord3d gradient (eval_evec[0].second[3*j], eval_evec[0].second[3*j+1], eval_evec[0].second[3*j+2]);
				gradients.push_back(gradient / gradient.norm());
			}

			if (gradients.size() != coords.size()) {
				cerr << "gradient and coordinate vector are of different size" << endl;
				break;
			}

			vector<coord3d> displacement;

			threadstream << "negative deflection" << endl;
			for (vector<coord3d>::size_type	j = 0; j < coords.size(); j++) {
				displacement.push_back(coords[j] - (gradients[j] * 0.1));
			}

			KS.setCoordinates(displacement);

			newKS = KS.optimize(threadstream, algo_switch, potential_switch, p, opt);

			hessian = newKS.hessian(p);
			eigenValues = diag(hessian);
			newKS.setHessian (eigenValues);
			
			newKS.setNumber(notMinimumKS[i].getNumber());


			//Inertia
			coord3d CoM = newKS.centreOfMass();
			newKS.shiftToCoM(CoM);
			vector< vector<double> > inertiaTensor = newKS.momentOfInertia();
			vector<double> inertia = diag(inertiaTensor);
			newKS.setMomentOfInertia(inertia);

		}
		while (reopts < 6 && !newKS.isMinimum());
		
		if (!newKS.isMinimum()) {
			threadstream << "Structure did not converge after " << reopts << " attempts." << endl;
			cerr << "Reopt failure, structure " << newKS.getNumber() << ", negative displacement" << endl;
		}
		reoptKS.push_back(newKS);
		threadstream << "***************End of Opt***************" << endl;
		remin << threadstream.rdbuf() << endl;	
	}
	remin.close();
	cout << "\tTime for structure optimization: " << T.timing() << " s" << endl << endl;






	//
	//WRITE COORD CHECKPOINT
	
	stringstream out;
	simpleout(optKS, out);
	simpleout(reoptKS, out);

	ofstream coord;
	coord.precision(14);
	coord.open("coord");
	coord << out.rdbuf();
	coord.close();
	
	
	cout << "\tTotal rumtime: " << T.total_timing() << " s" << endl;
    return 0; 

    
}
