#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <libconfig.h++>
#include "structure.h"
#include "iop.h"
#include "parameter.h"
#include "globals.h"

using namespace std;


using namespace libconfig;

void xyzout (structure &outputStructure, const string &name) 
{
    vector<coord3d> coordinates = outputStructure.getCoordinates();
    ofstream xyz;
    xyz.open ("output/" + name);
    xyz << outputStructure.nAtoms() << endl;
    xyz << "Energy: " << outputStructure.getEnergy() << ", moment of inertia: " << outputStructure.getMomentOfInertia() << endl;
    for (int i = 0; i < outputStructure.nAtoms(); i++) {
        xyz << setw(5) << left << "X" << setw(15) << setprecision(6) <<
            setiosflags(ios::fixed) << right << coordinates[i][0] << setw(15)
            << coordinates[i][1] << setw(15) << coordinates[i][2] << endl;
    }
    xyz.close();
}


void xyzoutall (vector<structure> &outputStructures, const string &name)
{
    ofstream xyz;
    xyz.open("output/" + name);
    for (vector<structure>::size_type i = 0; i < outputStructures.size(); i++)
    {
        vector<coord3d> coordinates = outputStructures[i].getCoordinates();
        xyz << outputStructures[i].nAtoms() << endl;
        xyz << "Energy: " << outputStructures[i].getEnergy() 
            << ", moment of inertia: " << outputStructures[i].getMomentOfInertia() << endl;
        for (int j = 0; j < outputStructures[i].nAtoms(); j++) 
        {
            xyz << setw(5) << left << "X" << setw(15) << setprecision(6) <<
                setiosflags(ios::fixed) << right << coordinates[j][0] <<
                setw(15) << coordinates[j][1] << setw(15) << coordinates[j][2]
                << endl;
        }
    }
    xyz.close();
}



void simpleout (vector<structure> &outputStructures, stringstream &output) {
    output.precision(14);
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
        structure kissingSphere(number, coordinates); //kissingSphere refers to 1 cluster of kissing Spheres
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

/*-------------------------------------------------------------------------------------*/
//                        function to read settings file
//                        return value determins potential
//                        0 = fail
//                        1 = LJ
//                        2 = ELJ
/*-------------------------------------------------------------------------------------*/

int readsettings (parameter<double> &opt, vector<double> &p, parameter<int> &switches, double &scalingFactor) {

    libconfig::Config cfg;
    try {
        cfg.readFile("settings");
    }
    catch(const FileIOException &fioex) {
        cerr << "\tI/O error while reading file." << endl;
        return 0;
    }
    catch(const ParseException &pex) {
        cerr << "\tParse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << endl;
        return 0;
    }

    const Setting &root = cfg.getRoot();

    string potential, algo;
    double convergence(1e-3);
    int potential_switch(0), algo_switch, scaling_switch;
    int nsteps(100);


    cout << endl;

    if (cfg.lookupValue("scaling.factor", scalingFactor)) {
        cout << "\tScaling: " << scalingFactor << endl;
        scaling_switch = 1;
    }
    else {
        cout << "\tNo scaling" << endl;
    }

    if (cfg.lookupValue("potential.name", potential)) {
        if (potential == "LJ") 
        {
            potential_switch = 1;
        }
        else if (potential == "ELJ") 
        {
            potential_switch = 2;
        }
        else if (potential == "RangeLJ")
        {
            potential_switch = 3;
        }
        else
        {
            cerr << potential << " is not a valid name for a potential." << endl;
            return 0;
        }
        cout << "\tPotential: " << potential << endl;
    }
    else {
        cout << "\tNo 'potential' setting in configuration file." << endl;
        return 0;
    }

    
    if (potential == "LJ") {
        double exp1(12), exp2(6), rm, epsilon;
        if (cfg.lookupValue("potential.epsilon", epsilon)) {
            p.push_back(epsilon);
            cout << "\t\tEpsilon: " << epsilon << endl;
        }
        else {
            cout << "\tNo 'epsilon' in configuration file." << endl;
            return 0;
        }
        if (cfg.lookupValue("potential.rm", rm)) {
            p.push_back(rm);
            cout << "\t\tRm: " << rm << endl;
        }
        else {
            cout << "\tNo 'rm' in configuration file." << endl;
            return 0;
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
        return 0;
    }


    if (cfg.lookupValue("opt.convergence", convergence)) {
        cout << "\t\tConvergence:  " << convergence << endl;
    }
    else {
        cout << "\t\tConvergence:  " << convergence << endl;
    }


    if (cfg.lookupValue("opt.nsteps", nsteps)) {
        cout << "\t\tNsteps: " << nsteps << endl;
    }
    else {
        cout << "\t\tNsteps: " << nsteps << endl;
    }
    opt.set("convergence", convergence);   //opt[1] 
    opt.set("nsteps", nsteps);   //opt[3]

    switches.set("potential", potential_switch); 
    switches.set("algo", algo_switch);
    switches.set("scaling", scaling_switch);

    cout << endl;

    return potential_switch;
}




template <typename T> void matrixout (vector< vector<T> > &matrix, ostream &out)
{
    //ios::fmtflags f(out.flags());
    ios oldState(nullptr);
    oldState.copyfmt(out);
    out << fixed;
    for (typename vector< vector<T> >::size_type i = 0; i < matrix.size(); i++)
    {
        for (typename vector<T>::size_type j = 0; j < matrix[i].size(); j++)
        {
            if (j == matrix[i].size() - 1)
            {
                out << setw(3) << ( (matrix[i][j] == 0) ? setprecision(0) :
                        setprecision(4) ) << matrix[i][j] << "\\\\"; 
                continue;
            }
            out << setw(3) << ( (matrix[i][j] == 0) ? setprecision(0) :
                    setprecision(4) ) << matrix[i][j] << " &"; }
        out << endl;
    }
    out.copyfmt(oldState);
}

template void matrixout<int> (vector< vector<int> > &matrix, ostream &out);
template void matrixout<double> (vector< vector<double> > &matrix, ostream &out);

     
