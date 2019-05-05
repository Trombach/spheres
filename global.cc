#include <iostream>
#include "structure.h"
#include "basinhopping.h"
#include "acceptanceTest.h"
#include "timer.h"

using namespace std;

int main ()
{
    timer T;

    structure cluster(10); 

    BasinHopping hop(cluster);

    unique_ptr<AcceptanceTest> accept;
    Metropolis* metro = new Metropolis(10);
    accept.reset(metro);

    hop.run(accept);

    cout << "Time: " << T.total_timing() << " s" << endl;

    return 0;
}
