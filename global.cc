#include <iostream>
#include <cstdlib>
#include "structure.h"
#include "basinhopping.h"
#include "acceptanceTest.h"
#include "timer.h"

using namespace std;

int main (int argc, char *argv[])
{
    timer T;

    if (argc < 3)
    {
        cerr << "Please provide system size and number of steps. Terminating program." << endl;
        return 1;
    }

    int n = atoi(argv[1]);
    int steps = atoi(argv[2]);
    

    structure cluster(n); 

    vector<coord3d> coordinates = cluster.getCoordinates();
    //for (auto& i : coordinates) cout << i << endl;

    shared_ptr<AcceptanceTest> accept_shared;
    Metropolis* metro = new Metropolis();
    accept_shared.reset(metro);

    BasinHopping<StorageByEnergy> hop(cluster,accept_shared,steps);

    hop.run();

    cout << "N minima: " << hop.nStructures() << endl;
    hop.printEnergies(cout);


    cout << "Time: " << T.total_timing() << " s" << endl;

    return 0;
}
