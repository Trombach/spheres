#include <iostream>
#include "structure.h"
#include "basinhopping.h"
#include "acceptanceTest.h"
#include "timer.h"

using namespace std;

int main ()
{
    timer T;

    structure cluster(8); 


    shared_ptr<AcceptanceTest> accept_shared;
    Metropolis* metro = new Metropolis(10);
    accept_shared.reset(metro);

    BasinHopping hop(cluster,accept_shared,1000);

    hop.run();

    cout << hop.nStructures() << endl;


    cout << "Time: " << T.total_timing() << " s" << endl;

    return 0;
}
