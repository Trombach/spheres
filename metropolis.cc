#include <random>
#include <math.h>
#include <iostream>
#include "acceptanceTest.h"

using namespace std;

bool Metropolis::operator() (double oldEnergy, double newEnergy)
{
    if (newEnergy < oldEnergy) return true;
    
    random_device r;
    default_random_engine generator{r()};
    
    double beta = (newEnergy - oldEnergy) / _temp; 
    double w = exp(-beta);
    
    uniform_real_distribution<double> distribution(0,1.0);
    double random = distribution(generator);

    
    if (random > w) return false;

    return true;
}
