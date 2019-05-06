#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include "basinhopping.h"
#include "geometry.h"
#include "potential.h"
#include "iop.h"
#include "parameter.h"
#include "lina.h"

using namespace std;


int BasinHopping::run () 
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
        dummy.open("out");
        _currentStep = potential->optimize(dummy, _currentStep, switches, opt); 
        dummy.close();
        
        vector< vector<double> >hessian = potential->calcHessian(_currentStep);
        vector<double> eigenValues = diag(hessian);
        _currentStep.setHessian(eigenValues);

        if (_currentStep.isMinimum()) //check for true minimum
        {
            if (this->checkConf()) //check spherical container
            {
                double oldE, newE = _currentStep.getEnergy();

                if (_iteration == 1) {oldE = newE;}
                else {oldE = _previousStep.getEnergy();}

                if (this->acceptStep(oldE, newE)) //accept the step
                {
                    accepted = true;
                    _uniqueStructures.addCluster(_currentStep);
                }
                else {_currentStep = _previousStep;} //rejected
            }
            else{_currentStep = _previousStep;} //rejected because spherical container
        }


        cout.precision(5);
        cout << fixed << left << setprecision(10);
        cout << setfill(' ')  << left <<
            _iteration << left <<
            " oldE = " << left << _previousStep.getEnergy() <<  left <<
            " newE = " << left << _currentStep.getEnergy() << left <<
            " accepted = " << boolalpha << accepted <<
            endl;

        this->propagate();

    }

    int n(0);
    for (auto& i : _uniqueStructures)
    {
        cout << i.first << endl;
        xyzout(i.second, to_string(n) + ".xyz");
        n++;
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

bool BasinHopping::acceptStep (double oldE, double newE)
{
    return (*_accept)(oldE, newE);
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
