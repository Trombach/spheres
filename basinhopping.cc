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
#include "storage.h"

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
    
    ofstream dummy;
    //dummy.open("out");

    float progress;

    for (int i = 0; i < _nsteps; i++)
    {
        _accepted = false;
        _iteration++;
        _previousStep.setEnergy(potential->calcEnergy(_previousStep)); 
        _currentStep.setEnergy(potential->calcEnergy(_currentStep)); 

        _currentStep = potential->optimize(dummy, _currentStep); 

        if (!_currentStep.isConverged()) 
        {
            _currentStep = _previousStep; //rejected because optimisation not converged
        }
        else
        {
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
                        _accepted = true;
                        _uniqueStructures.addCluster(_currentStep);
                    }
                    else {_currentStep = _previousStep;} //rejected
                }
                else{_currentStep = _previousStep;} //rejected because spherical container
            }
            else {_currentStep = _previousStep;} //rejected because of hessian
        }


        //cout.precision(5);
        //cout << fixed << left << setprecision(10);
        //cout << setfill(' ')  << left <<
        //    _iteration << left <<
        //    " oldE = " << left << _previousStep.getEnergy() <<  left <<
        //    " newE = " << left << _currentStep.getEnergy() << left <<
        //    " accepted = " << boolalpha << accepted <<
        //    endl;
        

        int barWidth = 70;
        progress = i / static_cast<double>(_nsteps);
        cout << "[";
        int pos = barWidth * progress;
        for (int j = 0; j < barWidth; j++)
        {
            if (j < pos) cout << "=";
            else if (j == pos) cout << ">";
            else cout << " ";
        }
        cout << "] " << static_cast<int>(progress * 100.0) << " % N: " << this->nStructures() << " step: " << _stepScale << " T: " << _accept->getT() << "\r";
        cout.flush();

        this->updateStep();
        this->propagate();


    }
    cout << endl;

    //dummy.close();

    _uniqueStructures.printStructures();


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
            i[j] += _stepScale * distribution(generator);
        }
    }

    _currentStep.setCoordinates(coordinates);
}

void BasinHopping::updateStep()
{
    if (_iteration == 1) return; 

    _nattempts++;
    if (_accepted) _naccept++;

    if ( fabs(_currentStep.getEnergy() - _previousStep.getEnergy()) < 1e-4) _nsame++;

    if (_nattempts % _interval == 0)
    {
        adjustStep();
        adjustTemp();
        resetUpdateStep(); 
    }
}

void BasinHopping::adjustStep ()
{
    double f = 1 - (static_cast<double>(_nsame) / _nattempts);
    if (f < 0.71) _stepScale /= 0.95; 
    else _stepScale *= 0.95; 
}


void BasinHopping::adjustTemp ()
{
    double temp = _accept->getT();
    int ndiff = _nattempts - _nsame;
    int ndiff_accept = _naccept - _nsame;

    double f(0);
    if (ndiff == 0) f = 1; 
    else f = static_cast<double>(ndiff_accept) / ndiff;

    if (f > 0.71) _accept->setT(temp * 0.95);
    else _accept->setT(temp / 0.95);
}



void BasinHopping::resetUpdateStep ()
{
    _nattempts = 0;
    _naccept = 0;
    _nsame = 0;
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
