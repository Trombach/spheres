#include <iostream>
#include <vector>
#include <map>
#include <functional>
#include <iomanip>
#include "../iop.h"
#include "../geometry.h"
#include "../structure.h"

using namespace std;

typedef map <double, unsigned int, function<bool(double a, double b)> > ncMap;

int main (int argc, char *argv[])
{

    cout << endl;

    string fileName = argv[argc - 1];

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

    cout << "\t+ Reading and processing structures" << endl;
    vector<structure> optKS = readallstruct(fileName);
    cout << "\t\t#structures :" << optKS.size() << endl;



/*-------------------------------------------------------------------------------------*/
/*                       output structures with nc = (40)                              */
/*-------------------------------------------------------------------------------------*/

    vector<structure> found;
    for (vector<structure>::size_type i = 0; i < optKS.size(); i++)
    {
        unsigned int nc(0);
        vector<double> interPartDist = optKS[i].getInterPartDist();
        for (vector< vector<double> >::size_type j = 0; j < interPartDist.size(); j++)
        {
            if (fabs(interPartDist[j] - 1) < 1e-8) 
            {
                nc++;
            }
        }
        //cout << "Structure " << optKS[i].getNumber() << " " << nc << endl; 
        
        if (nc == 36) found.push_back(optKS[i]);
       // matrixout(distMatrix, cout);
    }

    xyzoutall (found, "Nc.xyz");
    ofstream foundNc;
    stringstream out;
    foundNc.open("Nc");
    simpleout (found, out);
    foundNc.precision(14);
    foundNc << out.rdbuf();
    foundNc.close(); 



    return 0;
}
