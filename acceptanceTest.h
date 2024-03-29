#ifndef ACCEPTANCETEST
#define ACCEPTANCETEST

/*
Acceptance test for basinhopping step.  
*/

class AcceptanceTest
{
    protected:
        double _temp;

    public:
        AcceptanceTest(double temp) : _temp(temp) {}
        virtual ~AcceptanceTest() {};
        virtual bool operator() (double oldE, double newE) = 0;  
        void setT (double temp) {_temp = temp;}
        double getT() { return _temp;}

};

class Metropolis : public AcceptanceTest
{
    public:
        Metropolis(double temp = 5) : AcceptanceTest(temp) {}
        bool operator() (double oldEnergy, double newEnergy);
};

#endif
