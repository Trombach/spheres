#include <iostream>
#include "structure.h"
#include "iop.h"

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

    cout << endl;

    vector<structure> allKS = readallstruct(fileName);

    xyzoutall(allKS);

    return 0;
}

