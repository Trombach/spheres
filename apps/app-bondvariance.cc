#include <iostream>
#include <vector>
#include "../iop.h"
#include "../geometry.h"
#include "../structure.h"
#include "../potential.h"

using namespace std;

int main (int argc, char *argv[])
{
    cout << endl;


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

    for (vector<structure>::size_type i = 0; i < optKS.size(); i++)
    {
        double energy = potential->calcEnergy(optKS[i]);
        optKS[i].setEnergy(energy);
        vector<double> bonds;
        vector<coord3d> coords = optKS[i].getCoordinates();
        for (vector<coord3d>::size_type j = 0; j < coords.size(); j++)
        {
            for (vector<coord3d>::size_type k = j + 1; k < coords.size(); k++)
            {
                double dist = coord3d::dist(coords[j], coords[k]);
                double diff = fabs(dist - p[1]);
                if (diff < 0.1) 
                {
                    bonds.push_back(dist);
                }
            }
        }
    sort (bonds.begin(), bonds.end());
    cout << optKS[i].getNumber() << ": " << optKS[i].getEnergy() << " " << (bonds.back() - bonds.front()) << endl;
    }

        


    
    
   
    return 0;
}
