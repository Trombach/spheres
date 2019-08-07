#ifndef STRUCTURE
#define STRUCTURE

#include <vector>
#include "geometry.h"
#include "lina.h"
#include "parameter.h"
#include "graph.h"
#include "globals.h"

/*
This class handles 3d clusters made up of spheres.

Implementations can be found in the code file "structure.cc"

Parameters
----------

_energy :
    Energy of the cluster structure. Useful for sorting and analysis.
_number :
    A number that can be given to label the structure in a set of other
    structures.
_coordinates :
    This is the storage for the 3d coordinates. Each atoms coordinates are
    stored in an object of type "coord3d".
_momentOfInertia :
    Array that contains the moment of inertia eigenvalues.
_hessian :
    Array that contains the hessian eigenvalues.
_interPartDist :
    Array of all unique inter particle distances.
_adjMatrix :
    Holds connectivity information in terms of an adjacency matrix Aij = dij,
    with dij = 1 if atoms are said to be bound and 0 otherwise.
_distMatrix :
    Like _interPartDist but in matrix form (symmetric).
_bondVector :
    UNUSED! Was meant for analysis. 
_adjMatrix_eigenvalues :
    UNUSED! Eigenvalues of the adjacency matrix.
_uGraph :
    Undirected Graph object implemented with the boost graph library. Contains
    same information as adjacency matrix.
_converged :
    Boolean value that can be set to mark a structures as converged to a local
    minimum on the PES.

Constructors
------------
structure() :
    Default constructor.

structure ( int number, 
            std::vector<coord3d> coordinates, 
            bool calcProp = true) :
    Constructs instantiation with _number = number, _coordinates = coordinates
    and calculates properties if calcPror == true, which enables the member
    function "setCoordinates" to be used. The function shifts the cluster to
    its centre of mass, calculates the principal axes (by calculation the
    moment of inertia tensor and diagonalising it), rotates the cluster such
    that it is aligned with the principal axes and calculates the vector of
    inter particle distances (_interPartDist) and the distance matrix
    (_distMatrix).

structure (int number) :
    Initialises the cluster from random coordinates. Number of particles is
    equal to "number".
*/

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
    bool _converged;

public:   
    structure() :   _energy(0), 
                    _number(0), 
                    _coordinates(), 
                    _momentOfInertia {0,0,0},
                    _hessian(0),
                    _interPartDist(),
                    _bondVector(),
                    _adjMatrix_eigenvalues(),
                    _converged(false)
    {}

    /* constructor calculates several properties based on coordinates on creation */    
    structure ( int number, 
                std::vector<coord3d> coordinates, 
                bool calcProp = true) : _energy(0),
                                        _number(number), 
                                        _hessian(0),
                                        _interPartDist(),
                                        _bondVector(),
                                        _adjMatrix_eigenvalues(),
                                        _converged(false)
    { 
        if (calcProp) this->setCoordinates(coordinates);
        else 
        {
            _coordinates = coordinates;
            _converged = false;
        }
    }

    structure ( int number) :   _energy(0),
                                _number(0),
                                _hessian(0),
                                _interPartDist(),
                                _bondVector(),
                                _adjMatrix_eigenvalues(),
                                _converged(false)
    {
        this->randomize(number);
    }    

    void randomize(int number);
    bool isConverged() {return _converged;}
    void setConverged() {_converged = true;}

    //getter functions    
    int getNumber() const { return _number; }
    double getEnergy() const { return _energy; }
    std::vector<coord3d> getCoordinates() const { return _coordinates; }
    column_vector getFlattenedCoordinates();
    std::vector<double> getMomentOfInertia() const { return _momentOfInertia; }
    std::vector<double> getHessian() const { return _hessian; }
    std::vector<double> getInterPartDist() const { return _interPartDist; }
    std::vector< std::vector<int> > getAdjMatrix() const { return _adjMatrix; }
    std::vector<int> getBondVector() const { return _bondVector; }
    std::vector<double> getAdjMatrix_eigenvalues() const { return _adjMatrix_eigenvalues; }
    std::vector< std::vector<double> > getDistMatrix() const { return _distMatrix; }
    undirectedGraph getGraph() const { return _uGraph; }

    //setter functions
    void setNumber (int number) { _number = number; }
    void setEnergy (double energy) { _energy = energy; }
    void setCoordinates (std::vector<coord3d> coordinates) 
    { 
        _coordinates = coordinates; 
        _converged = false;

        this->shiftToCoM();

        std::vector< std::vector<double> > inertiaTensor = this->momentOfInertia();
        _momentOfInertia = diag(inertiaTensor);

        matrix3d axis = this->m3d_principalAxis ();
        this->rotateToPrincipalAxis(axis);

        this->propertyInterPartDist();

        this->propertyDistMatrix();
    }
    void setHessian (std::vector<double> hessianEigenvalues) {_hessian = hessianEigenvalues; }

    //property* functions calculate a property and store it in the proper member variable
    void propertyInterPartDist();
    void propertyAdjMatrix (std::vector<double> &p);
    void propertyDistMatrix();
    void propertyGraph(double rm = 1, double eps = 1e-10);
    void propertyGraph_ignoreCenter(double rm = 1, double eps = 1e-10);

    int nAtoms() { return (this->getCoordinates()).size(); }
    static std::vector<coord3d> unflattenCoordinates (column_vector &v);

    //for scaling structures structure object can just be multiplied with a double
    structure &operator*= (const double &y) 
    {
        for (int i = 0; i < this->nAtoms(); i++) { (*this)[i] *= y; }
        return *this;
    }
    structure operator* (const double &y) const { return structure(*this) *= y; }
    //direct access to atomic coordinates via [] operator
    coord3d& operator[] (unsigned int i) { return _coordinates[i]; }
    bool operator< (const structure y) const { return (this->getEnergy() < y.getEnergy()); }
    bool compareCoordinates (structure &y) const;

    //Inertia and centre of mass calculations
    coord3d centreOfMass ();
    void shiftToCoM ();
    std::vector< std::vector<double> > momentOfInertia ();
    matrix3d m3d_momentOfInertia ();
    matrix3d m3d_principalAxis ();
    void rotateToPrincipalAxis (matrix3d &principalAxis);

    //functions that depend on hessian matrix
    bool isMinimum ();
    int countNegativeEval ();


    //graph stuff
    std::vector< std::vector<int> > createAdjMatrix (std::vector<double> &p);
    undirectedGraph createGraph (double rm = 1, double eps = 1e-10);
    undirectedGraph createGraph_ignoreCenter (double rm = 1, double eps = 1e-10);
    std::vector<int> createBondVector ();
    std::vector<double> createAdjMatrix_egenvalues ();

    //nearest neighbours
    double longest_nearest_neighbour_distance (unsigned int n);

    //symmetry operations
    std::vector<coord3d> sig (int a);
    std::vector<coord3d> c2 (int a, int b);
    std::vector<coord3d> inv ();

};

#endif
