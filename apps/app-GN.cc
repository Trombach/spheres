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
/*             output structures with 12 spheres touching center                       */
/*-------------------------------------------------------------------------------------*/

    auto compare_double = [&] (double a, double b)
    {
        double eps = 1e-3;
        if (b - a > eps) return true;
        return false;
    };
    
    ncMap contacts_map(compare_double);

    vector<structure> Nc12;
    vector<double> max;
    vector<int> Nc;

    for (vector<structure>::size_type i = 0; i < optKS.size(); i++)
    {
        vector< vector<double> > distMatrix = optKS[i].getDistMatrix();
        bool gregory(false);
        unsigned int totalnc(0);

        for (vector< vector<double> >::size_type j = 0; j < distMatrix.size(); j++)
        {
            unsigned int nc(0);
            for (vector<double>::size_type k = 0; k < distMatrix[j].size(); k++)
            {
                if (fabs(distMatrix[j][k] - 1) < 1e-10) 
                {
                    nc++;
                    totalnc++;
                }
            }
            if (nc == 12)
            {
                gregory = true;
                double max_dist = *max_element(distMatrix[j].begin(), distMatrix[j].end());

                max.push_back(max_dist);
                if (max_dist < 1.376)
                {
                    Nc12.push_back(optKS[i]);
                }

                map<double, unsigned int, function<bool(double a, double b)> >::iterator iter(contacts_map.find(max_dist));
                if (iter != contacts_map.end())
                {
                    iter->second++;
                }
                else
                {
                    contacts_map[max_dist] = 1;
                }

                //matrixout(distMatrix, cout);
            }
        }
        if (gregory) Nc.push_back(totalnc / 2);
    }


    //for (map<double, unsigned int, function<bool(double a, double b)> >::iterator iter = contacts_map.begin(); iter != contacts_map.end(); iter++) 
    //{
    //    cout << setw(10) << right << iter->first << setw(30) << iter->second << setw(10) << endl;
    //}
    
    
    if (Nc.size() != max.size())
    {
        cerr << "Something went wrong" << endl;
        return 1;
    }
    for (vector<double>::size_type i = 0; i < max.size(); i++)
    {
        cout << max[i] << " " << Nc[i] << endl;
    }

    cout << "min of max " << setprecision(20) << *min_element(max.begin(), max.end()) << endl;
    xyzoutall (Nc12, "Nc12.xyz");
    ofstream foundNc12;
    stringstream out1;
    foundNc12.open("Nc12");
    simpleout (Nc12, out1);
    foundNc12.precision(14);
    foundNc12 << out1.rdbuf();
    foundNc12.close(); 


    return 0;
}
