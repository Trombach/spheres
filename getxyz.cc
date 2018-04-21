#include <iostream>
#include "structure.h"
#include "iop.h"
#include "potential.h"

using namespace std;

int main (int argc, char *argv[])
{
    cout << endl;

    if (argc < 2)
    {
        cerr << "Please provide filename!" << endl;
        return 1;
    }


    string fileName = argv[argc-1];

    if (fexists(fileName))
    {
        cout << "\tFile " << fileName << " exists." << endl;
    }
    else
    {
        cout << "\tFile " << fileName << " does not exist." << endl;
        return 1;
    }

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

    vector<structure> allKS = readallstruct(fileName);

    for (vector<structure>::size_type i = 0; i < allKS.size(); i++)
    {
        allKS[i].setEnergy(potential->calcEnergy(allKS[i]));
    }


    xyzoutall(allKS);

    return 0;
}

