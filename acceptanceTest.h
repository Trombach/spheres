#ifndef ACCEPTANCETEST
#define ACCEPTANCETEST

class AcceptanceTest
{
    public:
        virtual ~AcceptanceTest() {};
        virtual bool operator() (double oldE, double newE) = 0;  
};

class Metropolis : public AcceptanceTest
{
    private:
        double _temp;

    public:
        Metropolis(double temp) : _temp(temp) {}
        bool operator() (double oldEnergy, double newEnergy);
};

#endif
