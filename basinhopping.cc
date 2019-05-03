#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include "basinhopping.h"
#include "geometry.h"
#include "potential.h"
#include "iop.h"
#include "parameter.h"

using namespace std;


int BasinHopping::run (unique_ptr<AcceptanceTest>& accept) 
{

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
        case 3:
            potential.reset( RangeLJ::readPotential() );
            break;
        default:
            cerr << "readsettings returned an unknown status" << endl;
            return 1;
    }
    
    for (int i = 0; i < _nsteps; i++)
    {
        bool accepted(false);
        _iteration++;
        _previousStep.setEnergy(potential->calcEnergy(_previousStep)); 
        _currentStep.setEnergy(potential->calcEnergy(_currentStep)); 

        ofstream dummy;
        _currentStep = potential->optimize(dummy, _currentStep, switches, opt); 

        if (this->checkConf())
        {
            double oldE, newE = _currentStep.getEnergy();

            if (_iteration == 1) {oldE = newE;}
            else {oldE = _previousStep.getEnergy();}

            if (this->acceptStep(oldE, newE, accept)) 
            {
                accepted = true;
            }
            else
            {
                _currentStep = _previousStep;
            }
        }
        else
        {
            cerr << "Cluster not within spherical container" << endl;
            _currentStep = _previousStep;
        }

        cout << setfill(' ') << setw(10) << left <<
            _iteration << setw(10) << left <<
            _previousStep.getEnergy() << setw(10) << left <<
            _currentStep.getEnergy() << setw(10) << left <<
            boolalpha << accepted <<
            endl;

        this->propagate();

    }

    return 0;
}

void BasinHopping::propagate()
{
    _previousStep = _currentStep;
    vector<coord3d> coordinates = _currentStep.getCoordinates();

    random_device r;
    default_random_engine generator{r()};

    for (auto& i : coordinates)
    {
        for (int j = 0; j < 3; j++)
        {
            uniform_real_distribution<double> distribution(-1.0,1.0);
            i[j] += distribution(generator);
        }
    }

    _currentStep.setCoordinates(coordinates);
}

bool BasinHopping::acceptStep (double oldE, double newE, unique_ptr<AcceptanceTest>& accept)
{
    return (*accept)(oldE, newE);
}



bool BasinHopping::checkConf()
{   
    vector<coord3d> coordinates = _currentStep.getCoordinates();
    
    for (auto& i : coordinates)
    {
        if (i.norm() > 4) return false;
    }

    return true;
}
