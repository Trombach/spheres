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
#include "timer.h"
#include "parameter.h"
#include "potential.h"
#include "globals.h"


using namespace std; 


typedef map <pair < double, vector<double> >, unsigned int, function<bool( pair < double, vector<double> > a, pair < double, vector<double> > b)> > energyMap;


//MAIN FUNCTION BEGINS HERE

int main (int argc, char *argv[]) {
    timer T;



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
    parameter<double> opt; //vector of algo settings, 0 == accuracy, 1 == dforce, 2 == stepsize, 3 == nsteps
    parameter<int> switches; //vector of switches, 0 == potential, 1 == algo, 2 == scaling
    double scalingFactor(1.0);

    int status = readsettings(opt, p, switches, scalingFactor);
    std::unique_ptr< pairPotential > potential;
    switch (status)
    {
        case 0:
            cerr << "Failed reading settings file" << endl;
            return 1;
        case 1:
            potential.reset( LJ::readPotential() );
            break;
        case 2:
            potential.reset( ELJ::readPotential() );
            break;
        case 3:
            potential.reset( RangeLJ::readPotential() );
            break;
        default:
            cerr << "readsettings returned an unknown status" << endl;
            return 1;
    }
    
    
    cout << endl;






    
    //READ IN ALL STRUCTURES AT ONCE
    cout << "\t+ Reading and processing structures" << endl;
    vector<structure> optKS = readallstruct(fileName);
    cout << "\t\t#structures :" << optKS.size() << endl;
    
    //if scaling is found in settings file scale all coordinates accordingly
    const int scaling_switch = switches.get("scaling");
    switch (scaling_switch) {
        case 1:
            for (vector<structure>::size_type i = 0; i < optKS.size(); i++) {
                optKS[i] *= scalingFactor;
            }
        default:
            break;
    }






    //HESSIAN, etc
    ofstream energystats;
    unsigned int hessianWarnings = 0;
    vector<int> notMinimum;
    vector<structure> notMinimumKS;

    #pragma omp parallel for
    for (vector<structure>::size_type i = 0; i < optKS.size(); i++) {

        vector< vector<double> > hessian;
        vector<double> eigenValues;
        stringstream threadstream;

        optKS[i].setEnergy(potential->calcEnergy(optKS[i]));

        hessian = potential->calcHessian(optKS[i]);
        eigenValues = diag(hessian);
        optKS[i].setHessian (eigenValues);

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

    //remove non-min structures from optKS
    for (unsigned int i = 0; i < notMinimum.size(); i++)
    {
        optKS.erase(remove_if(optKS.begin(), optKS.end(), [&] (const structure &s)
                -> bool {return (s.getNumber() == notMinimum[i]);}), optKS.end());
    }


    cout << "\t\t#non-minimum structures: " << notMinimumKS.size() << endl;
    cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;



    vector<structure> optKS0;
    optKS0.push_back(optKS[0]);




    
    //compare function for vectors of doubles
    auto compare_vector_double = [&] (vector<double> a, vector<double> b) {
        assert (a.size() == b.size());
        double eps = 1e-5;
        vector<double> diff;
        for (vector<double>::size_type i = 0; i < a.size(); i++) {
            diff.push_back(abs(a[i]-b[i]));
        }
        auto max = max_element (begin(diff), end(diff));
        if (*max < eps) return true;
        return false;
    };









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
        double eps = 1e-4;
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

    energystats << endl << "Number of warnings for Hessian eigenvalues: " << hessianWarnings << endl;
    energystats << "List of non-minimum Structures: " << notMinimum << endl;

    energies.close();
    energystats.close();



    
    cout << "\t+ Energy and moment of inertia analysis:" << endl;
    cout << "\t\tUnique structures: " << energyStat.size() << endl;
    cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;











    //sort all minimum structures into vector< vector <structure> > according to distance vector
    vector< vector<structure> > eqClasses_dist;
    eqClasses_dist.push_back(optKS0);
    for (vector<structure>::size_type i = 1; i < optKS.size(); i++) {
        bool matched(0);
        for (vector< vector<structure> >::size_type j = 0; j < eqClasses_dist.size(); j++) {
            if (compare_vector_double (optKS[i].getInterPartDist(), eqClasses_dist[j][0].getInterPartDist())) {
                eqClasses_dist[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
        }
        if (matched == 0) {
            vector<structure> newEqClass;
            newEqClass.push_back(optKS[i]);
            eqClasses_dist.push_back(newEqClass);
        }
    }


    cout << "\t+ Particle distance analysis" << endl;
    cout << "\t\tEquality classes: " << eqClasses_dist.size() << endl;
    cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;

    //output of groups for matching
    vector<structure> outputStructures;
    for (vector<structure> i : eqClasses_dist)
    {
        outputStructures.push_back(i[0]);
    } 

    stringstream out;
    simpleout(outputStructures, out); 

    ofstream optCoord;
    optCoord.precision(14);
    optCoord.open("optCoord");
    optCoord << out.rdbuf();
    optCoord.close();


    
    
    cout << "\t+ Total rumtime: " << T.total_timing() << " s" << endl;
    return 0; 

    
}

//END













//OLD FUNCTIONS NOT IN USE

/*
    //SORT BY direct structure comparison
    vector< vector<structure> > eqClasses_coord;
    eqClasses_coord.push_back(optKS0);
    for (vector<structure>::size_type i = 1; i < optKS.size(); i++) {
        bool matched(0);
        for (vector< vector<structure> >::size_type j = 0; j < eqClasses_coord.size(); j++) {
            if (optKS[i].compareCoordinates(eqClasses_coord[j][0])) {
                eqClasses_coord[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
            structure sig0;
            vector<coord3d> coord_sig0 = optKS[i].sig(0);
            sig0.setCoordinates(coord_sig0);
            if (sig0.compareCoordinates(eqClasses_coord[j][0])) {
                eqClasses_coord[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
            structure sig1;
            vector<coord3d> coord_sig1 = optKS[i].sig(1);
            sig1.setCoordinates(coord_sig1);
            if (sig1.compareCoordinates(eqClasses_coord[j][0])) {
                eqClasses_coord[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
            structure sig2;
            vector<coord3d> coord_sig2 = optKS[i].sig(2);
            sig2.setCoordinates(coord_sig2);
            if (sig2.compareCoordinates(eqClasses_coord[j][0])) {
                eqClasses_coord[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
            structure c201;
            vector<coord3d> coord_c201 = optKS[i].c2(0,1);
            c201.setCoordinates(coord_c201);
            if (c201.compareCoordinates(eqClasses_coord[j][0])) {
                eqClasses_coord[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
            structure c202;
            vector<coord3d> coord_c202 = optKS[i].c2(0,2);
            c202.setCoordinates(coord_c202);
            if (c202.compareCoordinates(eqClasses_coord[j][0])) {
                eqClasses_coord[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
            structure c212;
            vector<coord3d> coord_c212 = optKS[i].c2(1,2);
            c212.setCoordinates(coord_c212);
            if (c212.compareCoordinates(eqClasses_coord[j][0])) {
                eqClasses_coord[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
            structure inverse;
            vector<coord3d> coord_inv = optKS[i].inv();
            inverse.setCoordinates(coord_inv);
            if (inverse.compareCoordinates(eqClasses_coord[j][0])) {
                eqClasses_coord[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
            
        }
        if (matched == 0) {
            vector<structure> newEqClass;
            newEqClass.push_back(optKS[i]);
            eqClasses_coord.push_back(newEqClass);
        }
    }



    cout << "\t+ Comparison of coordinates" << endl;
    cout << "\t\tEquality classes: " << eqClasses_coord.size() << endl;
    cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;
*/


/*
    vector< vector<structure> > eqClasses_graph;
    eqClasses_graph.push_back(optKS0);
    for (vector<structure>::size_type i = 1; i < optKS.size(); i++) {
        bool matched(0);
        structure::undirectedGraph g0 = optKS[i].createGraph(p);
        for (vector<structure>::size_type j = 0; j < eqClasses_graph.size(); j++) {
            structure::undirectedGraph g1 = optKS[j].createGraph(p);

            typename boost::property_map<structure::undirectedGraph, boost::vertex_index_t>::type
            v0_index_map = boost::get (boost::vertex_index, g0);
            //v1_index_map = boost::get (boost::vertex_index, test1);

            vector<boost::graph_traits<structure::undirectedGraph>::vertex_descriptor> f(num_vertices(g0));
    
            if (boost::isomorphism (g0, g1, boost::isomorphism_map (make_iterator_property_map(f.begin(), v0_index_map, f[0])))) {
                eqClasses_graph[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
        }
        if (matched == 0) {
            vector<structure> newEqClass;
            newEqClass.push_back(optKS[i]);
            eqClasses_graph.push_back(newEqClass);
        }

    }


    cout << "\t+ Graph isomorphism" << endl;
    cout << "\t\tEquality classes: " << eqClasses_graph.size() << endl;
    cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;
*/


/*
    //post analysis of distance vector
    vector< vector<structure> > eqClasses_coord_dist;
    eqClasses_coord_dist.push_back(eqClasses_coord[0]);
    for (vector< vector<structure> >::size_type i = 1; i < eqClasses_coord.size(); i++) {
        bool matched(0);
        for (vector< vector<structure> >::size_type j = 0; j < eqClasses_coord_dist.size(); j++) {
            if (compare_vector_double (eqClasses_coord[i][0].getInterPartDist(), eqClasses_coord_dist[j][0].getInterPartDist())) {
                eqClasses_coord_dist[j].insert (eqClasses_coord_dist[j].end(), eqClasses_coord[i].begin(), eqClasses_coord[i].end());
                matched = 1;
                break;
            }
        }
        if (matched == 0) {
            eqClasses_coord_dist.push_back(eqClasses_coord[i]);
        }
    }


    vector<structure> structure0;
    for (vector< vector<structure> >::size_type i = 0; i < eqClasses_coord_dist.size(); i++) {
        structure0.push_back(eqClasses_coord_dist[i][0]);
        vector<double> inertia = eqClasses_coord_dist[i][0].getMomentOfInertia();
        stringstream in;
        in << inertia << endl;
        xyzout (eqClasses_coord[i][0], in.str() + to_string(i) + ".xyz");
    }

    cout << "\t+ Comparison of coordinates and distance vector" << endl;
    cout << "\t\tEquality classes: " << eqClasses_coord_dist.size() << endl;
    cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;
*/  

/*
    //sort all minimum structures into vector< vector<structure> > according to squared adj matrix
    vector< vector<structure> > eqClasses_adjMat2;
    eqClasses_adjMat2.push_back(optKS0);
    for (vector<structure>::size_type i = 1; i < optKS.size(); i++) {
        bool matched(0);
        for (vector< vector<structure> >::size_type j = 0; j < eqClasses_adjMat2.size(); j++) {
            if (optKS[i].getBondVector() == eqClasses_adjMat2[j][0].getBondVector()) {
                eqClasses_adjMat2[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
        }
        if (matched == 0) {
            vector<structure> newEqClass;
            newEqClass.push_back(optKS[i]);
            eqClasses_adjMat2.push_back(newEqClass);
        }
    }


    cout << "\t+ Bond number analysis" << endl;
    cout << "\t\tEquality classes: " << eqClasses_adjMat2.size() << endl;
    cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;
*/

/*  
    //sort all minimum structures into vector< vector<structure> > according to adj matrix eigenvalues
    vector< vector<structure> > eqClasses_adjMatEv;
    eqClasses_adjMatEv.push_back(optKS0);
    for (vector<structure>::size_type i = 1; i < optKS.size(); i++) {
        bool matched(0);
        for (vector< vector<structure> >::size_type j = 0; j < eqClasses_adjMatEv.size(); j++) {
            if (compare_vector_double (optKS[i].getAdjMatrix_eigenvalues(), eqClasses_adjMatEv[j][0].getAdjMatrix_eigenvalues())) {
                eqClasses_adjMatEv[j].push_back(optKS[i]);
                matched = 1;
                break;
            }
        }
        if (matched == 0) {
            vector<structure> newEqClass;
            newEqClass.push_back(optKS[i]);
            eqClasses_adjMatEv.push_back(newEqClass);
        }
    }


    cout << "\t+ Adjacency matrix eigen analysis" << endl;
    cout << "\t\tEquality classes: " << eqClasses_adjMatEv.size() << endl;
    cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;
*/
