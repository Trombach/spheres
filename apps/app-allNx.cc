#include <iostream>
#include <vector>
#include <map>
#include <functional>
#include <iomanip>
#include "../iop.h"
#include "../geometry.h"
#include "../structure.h"

using namespace std;

struct NxCollection
{
    unsigned int _totalnc, _nc3, _nc4, _nc5, _nc6, _nc7, _nc12;

    NxCollection() :    _totalnc(0),
                        _nc3(0),
                        _nc4(0),
                        _nc5(0),
                        _nc6(0),
                        _nc7(0),
                        _nc12(0)
    {}
};

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
/*                       find kissing number for each sphere                           */
/*-------------------------------------------------------------------------------------*/

    vector<NxCollection> allNx;
    for (vector<structure>::size_type i = 0; i < optKS.size(); i++)
    {
        unsigned int totalnc(0);
        unsigned int nc3(0), nc4(0), nc5(0), nc6(0), nc7(0), nc12(0);
        NxCollection Nx;
        vector< vector<double> > distMatrix = optKS[i].getDistMatrix();
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
            if (nc == 12) nc12++;
            if (nc == 3) nc3++;
            if (nc == 4) nc4++;
            if (nc == 5) nc5++;
            if (nc == 6) nc6++;
            if (nc == 7) nc7++;
        }
        Nx._totalnc = totalnc / 2;
        Nx._nc7 = nc7;
        Nx._nc6 = nc6;
        Nx._nc5 = nc5;
        Nx._nc4 = nc4;
        Nx._nc3 = nc3;
        Nx._nc12 = nc12;
        allNx.push_back(Nx);
    }

    auto ncCompare = [&] (NxCollection a, NxCollection b)
    {
        return  (a._totalnc > b._totalnc) ||
                ((a._totalnc == b._totalnc) && (a._nc6 > b._nc6)) ||
                ((a._totalnc == b._totalnc) && (a._nc6 == b._nc6) && (a._nc5 > b._nc6)) ||
                ((a._totalnc == b._totalnc) && (a._nc6 == b._nc6) && (a._nc5 == b._nc6) && (a._nc4 > b._nc4)) ||
                ((a._totalnc == b._totalnc) && (a._nc6 == b._nc6) && (a._nc5 == b._nc6) && (a._nc4 == b._nc4) && (a._nc3 > b._nc3));
    };

    sort(allNx.begin(), allNx.end(), ncCompare);

    //for (vector<NxCollection>::size_type i = 0; i < allNx.size(); i++)
    //{
    //    cout    << allNx[i]._totalnc << " "
    //            << allNx[i]._nc6 << " "
    //            << allNx[i]._nc5 << " "
    //            << allNx[i]._nc4 << " "
    //            << allNx[i]._nc3
    //            << endl;
    //}

    auto NxCollectionCompare = [&] (NxCollection a, NxCollection b)
    {
        return (a._totalnc == b._totalnc) && (a._nc6 == b._nc6) && (a._nc5 == b._nc5) && (a._nc4 == b._nc4) && (a._nc3 == b._nc3);
    };

    vector< vector<NxCollection> > mapping;
    vector<NxCollection> allNx0;
    allNx0.push_back(allNx[0]);
    mapping.push_back(allNx0);
    for (vector<NxCollection>::size_type i = 1; i < allNx.size(); i++)
    {
        bool matched(false);
        for (vector< vector<NxCollection> >::size_type j = 0; j < mapping.size(); j++)
        {
            if (NxCollectionCompare (allNx[i], mapping[j][0]))
            {
                mapping[j].push_back(allNx[i]);
                matched = true;
                break;
            }
        }
        if (matched == false)
        {
            vector<NxCollection> newEqClass;
            newEqClass.push_back(allNx[i]);
            mapping.push_back(newEqClass);
        }
    }
    
    for (vector< vector<NxCollection> >::size_type i = 0; i < mapping.size(); i++)
    {
        cout << mapping[i][0]._totalnc << " " << mapping[i][0]._nc7 << " " << mapping[i][0]._nc6 << " " << mapping[i][0]._nc5 << " " << mapping[i][0]._nc4 << " " << mapping[i][0]._nc3 << "   " << mapping[i].size() << endl;
    }





    return 0;
}
