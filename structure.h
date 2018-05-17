#ifndef STRUCTURE
#define STRUCTURE

#include <vector>
#include "geometry.h"
#include "lina.h"
#include "parameter.h"
#include "graph.h"


class structure {    
    
private:
    double _energy;
    int _number;
    std::vector<coord3d> _coordinates;
    std::vector<double> _momentOfInertia;
    std::vector<double> _hessian;
    std::vector<double> _interPartDist;
    std::vector< std::vector<int> > _adjMatrix;
    std::vector< std::vector<double> > _distMatrix;
    std::vector<int> _bondVector;
    std::vector<double> _adjMatrix_eigenvalues;
    undirectedGraph _uGraph;

public:   
    structure() :   _energy(0), 
                    _number(0), 
                    _coordinates(), 
                    _momentOfInertia {0,0,0},
                    _hessian(0),
                    _interPartDist(),
                    _bondVector(),
                    _adjMatrix_eigenvalues()
    {}

    /* constructor calculates several properties based on coordinates on creation */    
    structure (int number, std::vector<coord3d> coordinates) :  _energy(0),
                                                                _number(number), 
                                                                _coordinates(coordinates), 
                                                                _hessian(0),
                                                                _interPartDist(),
                                                                _bondVector(),
                                                                _adjMatrix_eigenvalues()
    { 
        this->shiftToCoM();

        std::vector< std::vector<double> > inertiaTensor = this->momentOfInertia();
        _momentOfInertia = diag(inertiaTensor);

        matrix3d axis = this->m3d_principalAxis ();
        this->rotateToPrincipalAxis(axis);

        this->propertyInterPartDist();

        this->propertyDistMatrix();
    }

    

    int getNumber() const { return _number; }
    double getEnergy() const { return _energy; }
    std::vector<coord3d> getCoordinates() const { return _coordinates; }
    std::vector<double> getMomentOfInertia() const { return _momentOfInertia; }
    std::vector<double> getHessian() const { return _hessian; }
    std::vector<double> getInterPartDist() const { return _interPartDist; }
    std::vector< std::vector<int> > getAdjMatrix() const { return _adjMatrix; }
    std::vector<int> getBondVector() const { return _bondVector; }
    std::vector<double> getAdjMatrix_eigenvalues() const { return _adjMatrix_eigenvalues; }
    std::vector< std::vector<double> > getDistMatrix() const { return _distMatrix; }
    undirectedGraph getGraph() const { return _uGraph; }

    void setNumber (int number) { _number = number; }
    void setEnergy (double energy) { _energy = energy; }
    void setCoordinates (std::vector<coord3d> coordinates) 
    { 
        _coordinates = coordinates; 
        this->shiftToCoM();

        std::vector< std::vector<double> > inertiaTensor = this->momentOfInertia();
        _momentOfInertia = diag(inertiaTensor);

        matrix3d axis = this->m3d_principalAxis ();
        this->rotateToPrincipalAxis(axis);

        this->propertyInterPartDist();

        this->propertyDistMatrix();
    }
    void setMomentOfInertia (std::vector<double> inertiaEigenvalues) { _momentOfInertia = inertiaEigenvalues; }
    void setHessian (std::vector<double> hessianEigenvalues) {_hessian = hessianEigenvalues; }
    void setInterPartDist (std::vector<double> distances) { _interPartDist = distances;}
    void setAdjMatrix (std::vector< std::vector<int> > adjMatrix) { _adjMatrix = adjMatrix; }
    void setBondVector (std::vector<int> bondVector) { _bondVector = bondVector; }
    void setAdjMatrix_eigenvalues (std::vector<double> eigenvalues) { _adjMatrix_eigenvalues = eigenvalues; }
    void setDistMatrix (std::vector< std::vector<double> > distMatrix) { _distMatrix = distMatrix; }
    void setGraph (undirectedGraph uGraph) { _uGraph = uGraph; }


    void propertyInterPartDist();
    void propertyAdjMatrix (std::vector<double> &p);
    void propertyDistMatrix();
    void propertyGraph(double rm = 1, double eps = 1e-10);
    void propertyGraph_ignoreCenter(double rm = 1, double eps = 1e-10);

    int nAtoms() { return (this->getCoordinates()).size(); }

    void push_back(coord3d spheres) { _coordinates.push_back(spheres); }

    structure &operator*= (const double &y) 
    {
        for (int i = 0; i < this->nAtoms(); i++) { (*this)[i] *= y; }
        return *this;
    }
    structure operator* (const double &y) const { return structure(*this) *= y; }
    coord3d& operator[] (unsigned int i) { return _coordinates[i]; }
    bool operator< (const structure y) const { return (this->getEnergy() < y.getEnergy()); }
    bool compareCoordinates (structure &y) const;


/*--------------------------------------------------------------------------------------*/
/*                           old optimizer stuff                                        */
/*--------------------------------------------------------------------------------------*/

    //function to sum over all sphere interactions, change later to work with different potentials
    double sumOverAllInteractions (const std::vector<double> &p);
    //function to sum over all gradients to get gradient for each sphere
    std::vector<coord3d> sumOverAllGradients (const std::vector<double> &p);
    //initialize minimizer function
    structure optimize (std::ostream &min, parameter<int> switches, std::vector<double> parameters, parameter<double> opt);

    std::vector< std::vector<double> > hessian (const std::vector<double> &p);

/*--------------------------------------------------------------------------------------*/


    coord3d centreOfMass ();
    void shiftToCoM ();
    std::vector< std::vector<double> > momentOfInertia ();
    matrix3d m3d_momentOfInertia ();
    bool isMinimum ();
    int countNegativeEval ();

    matrix3d m3d_principalAxis ();
    void rotateToPrincipalAxis (matrix3d &principalAxis);

    //graph stuff
    std::vector< std::vector<int> > createAdjMatrix (std::vector<double> &p);
    undirectedGraph createGraph (double rm = 1, double eps = 1e-10);
    undirectedGraph createGraph_ignoreCenter (double rm = 1, double eps = 1e-10);
    std::vector<int> createBondVector ();
    std::vector<double> createAdjMatrix_egenvalues ();

    //nearest neighbours
    template <unsigned int n> double longest_nearest_neighbour_distance ();

    //symmetry
    std::vector<coord3d> sig (int a);
    std::vector<coord3d> c2 (int a, int b);
    std::vector<coord3d> inv ();

};

#endif
