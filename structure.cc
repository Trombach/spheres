#include <random>
#include "structure.h"
#include "lina.h"
#include "graph.h"
#include "globals.h"

using namespace std;
using namespace boost;

coord3d structure::centreOfMass () {
    double totalMass = nAtoms(); 
    coord3d coordinateSum;
    for (int i = 0; i < nAtoms(); i++) {
        coordinateSum += (*this)[i];
    }
    return coordinateSum * (1 / totalMass);
}


void structure::shiftToCoM () {
    coord3d CoM = this->centreOfMass();
    vector<coord3d> shiftedCoordinates;
    for (int i = 0; i < nAtoms(); i++) {
        shiftedCoordinates.push_back((*this)[i] - CoM);
    }
    _coordinates = shiftedCoordinates;
}

vector< vector<double> > structure::momentOfInertia () {
    vector< vector<double> > inertiaMatrix (3, vector<double> (3, 0));
    for (int i = 0; i < nAtoms(); i++) {
        inertiaMatrix[0][0] += pow((*this)[i][1],2) + pow((*this)[i][2],2);
        inertiaMatrix[1][1] += pow((*this)[i][2],2) + pow((*this)[i][0],2);
        inertiaMatrix[2][2] += pow((*this)[i][0],2) + pow((*this)[i][1],2);
        inertiaMatrix[0][1] += - (*this)[i][0] * (*this)[i][1];
        inertiaMatrix[0][2] += - (*this)[i][0] * (*this)[i][2];
        inertiaMatrix[1][2] += - (*this)[i][1] * (*this)[i][2];
        inertiaMatrix[1][0] += - (*this)[i][0] * (*this)[i][1];
        inertiaMatrix[2][0] += - (*this)[i][0] * (*this)[i][2];
        inertiaMatrix[2][1] += - (*this)[i][1] * (*this)[i][2];
    }
    return inertiaMatrix;
}

        
matrix3d structure::m3d_momentOfInertia () {
    matrix3d inertiaMatrix;
    for (int i = 0; i < nAtoms(); i++) {
        inertiaMatrix(0,0) += pow((*this)[i][1],2) + pow((*this)[i][2],2);
        inertiaMatrix(1,1) += pow((*this)[i][2],2) + pow((*this)[i][0],2);
        inertiaMatrix(2,2) += pow((*this)[i][0],2) + pow((*this)[i][1],2);
        inertiaMatrix(0,1) += - (*this)[i][0] * (*this)[i][1];
        inertiaMatrix(0,2) += - (*this)[i][0] * (*this)[i][2];
        inertiaMatrix(1,2) += - (*this)[i][1] * (*this)[i][2];
        inertiaMatrix(1,0) += - (*this)[i][0] * (*this)[i][1];
        inertiaMatrix(2,0) += - (*this)[i][0] * (*this)[i][2];
        inertiaMatrix(2,1) += - (*this)[i][1] * (*this)[i][2];
    }
    return inertiaMatrix;
}

matrix3d structure::m3d_principalAxis () {
    matrix3d I = this->m3d_momentOfInertia();
    matrix3d V = m3d_diagv (I);
    return V;
}

void structure::rotateToPrincipalAxis (matrix3d &principalAxis) {
    matrix3d transMatrix = principalAxis.inverse();
    //matrix3d test =  principalAxis * transMatrix;
    vector<coord3d> coord = this->getCoordinates(), newCoord;

    for (vector<coord3d>::size_type i = 0; i < coord.size(); i++) {
        newCoord.push_back (transMatrix * coord[i]);
    }

    _coordinates = newCoord;    
}   

column_vector structure::getFlattenedCoordinates()
{
    column_vector x((this->nAtoms()) * 3);
    for (int i = 0; i < this->nAtoms(); i++)
        for (int j = 0; j < 3; j++)
            x(3 * i + j) = (*this)[i][j];

    return x;
}

vector<coord3d> structure::unflattenCoordinates (column_vector &v)
{
    vector<coord3d> newCoordinates;
    for (int i = 0; i < v.size() / 3; i++)
    {
        coord3d sphere (v(3 * i), v(3 * i + 1), v(3 * i + 2));
        newCoordinates.push_back(sphere);
    }
    return newCoordinates;
}


bool structure::isMinimum () {
    vector<double> hessian = this->getHessian();
    if (hessian[0] > -0.001) return true;
    return false;
}

int structure::countNegativeEval ()
{
    int count(0);
    vector<double> hessian = this->getHessian();
    for (vector<double>::size_type i = 0; i < hessian.size(); i++)
    {
        if (hessian[i] < -0.001) count++;
    }
    return count;
}


undirectedGraph structure::createGraph (double rm, double eps) 
{
    vector< vector<double> > distMatrix = this->getDistMatrix();

    typedef pair<int, int> edge;
    vector<edge> edgeVec;
    for (vector< vector<double> >::size_type i = 0; i < distMatrix.size(); i++) 
    {
        for (vector< vector<double> >::size_type j = i + 1; j < distMatrix.size(); j++) 
        {
            double diff = fabs(distMatrix[i][j] - rm);
            if (diff < eps) edgeVec.push_back(edge(i,j));
        }
    }

    undirectedGraph g (edgeVec.begin(), edgeVec.end(), this->nAtoms());
    unsigned int i(0);
    BGL_FORALL_VERTICES_T(vertex, g, undirectedGraph)
    {
        g[vertex].coordinate = _coordinates[i];
        g[vertex].index = i;
        i++;
    }
    //cout << num_edges(g) << " " << num_vertices(g) << endl;

    return g;
}


undirectedGraph structure::createGraph_ignoreCenter (double rm, double eps) 
{
    vector< vector<double> > distMatrix = this->getDistMatrix();

    unsigned long centerIndex(0);
    for (vector< vector<double> >::size_type i = 0; i < distMatrix.size(); i++)
    {
        int nc(0);
        for (vector<double>::size_type j = 0; j < distMatrix[i].size(); j++)
        {
            if (fabs(distMatrix[i][j] - rm) < eps) nc++;
        }
        if (nc == 12) centerIndex = i;
    }

    typedef pair<int, int> edge;
    vector<edge> edgeVec;
    for (vector< vector<double> >::size_type i = 0; i < distMatrix.size(); i++) 
    {
        if (i == centerIndex) continue;
        for (vector< vector<double> >::size_type j = i + 1; j < distMatrix.size(); j++) 
        {
            if (j == centerIndex) continue;
            double diff = fabs(distMatrix[i][j] - rm);
            if (diff < eps) edgeVec.push_back(edge(i,j));
        }
    }

    undirectedGraph g (edgeVec.begin(), edgeVec.end(), this->nAtoms() - 1);
    unsigned int i(0);
    BGL_FORALL_VERTICES_T(vertex, g, undirectedGraph)
    {
        g[vertex].coordinate = _coordinates[i];
        g[vertex].index = i;
        i++;
    }
    //cout << num_edges(g) << " " << num_vertices(g) << endl;

    return g;
}

vector< vector<int> > structure::createAdjMatrix (vector<double> &p) {
    vector< vector<int> > adjMatrix;
    vector<coord3d> currentCoord = this->getCoordinates();
    for (vector<coord3d>::size_type j = 0; j < currentCoord.size(); j++) {
        vector<int> currentRow;
        for (vector<coord3d>::size_type k = 0; k < currentCoord.size(); k++) {
            if (j == k) currentRow.push_back(0);
            else {
                double diff = fabs(coord3d::dist (currentCoord[j], currentCoord[k]) - p[1]);
                if (diff < 0.1) currentRow.push_back(1);
                else currentRow.push_back(0);
            }
        }
        adjMatrix.push_back(currentRow);
    }
    return adjMatrix;
}




vector<int> structure::createBondVector () {
    vector< vector<int> > adjMatrix = this->getAdjMatrix();
    vector<int> bondVector (adjMatrix.size());
    for (vector< vector<int> >::size_type i = 0; i < adjMatrix.size(); i++) {
        for (vector< vector<int> >::size_type j = 0; j < adjMatrix.size(); j++) {
            bondVector[i] += adjMatrix[i][j] * adjMatrix[j][i];
        }
    }
    return bondVector;
}


vector<double> structure::createAdjMatrix_egenvalues () {
    vector< vector<int> > adjMatrix = this->getAdjMatrix();
    vector< vector<double> > adjMatrix_double;
    for (vector< vector<int> >::size_type i = 0; i < adjMatrix.size(); i++) {
        vector<double> row (adjMatrix[i].begin(), adjMatrix[i].end());
        adjMatrix_double.push_back(row);
    }
    vector<double> adjMatrix_eigenvalues = diag (adjMatrix_double);
    return adjMatrix_eigenvalues;
}

bool structure::compareCoordinates (structure &y) const {
    vector<coord3d> a = this->getCoordinates();
    vector<coord3d> b = y.getCoordinates();

    for (vector<coord3d>::size_type i = 0; i < a.size(); i++) {
        bool match(0);
        for (vector<coord3d>::size_type j = 0; j < b.size(); j++) {
            coord3d diff = a[i] - b[j];
            //cout << i << " " << diff.norm() << endl;
            if (diff.norm() < 1e-4) {
                match = 1;
                break;
            }
        }
        if (match == 0) return false;
    }
    return true;
}


void structure::propertyInterPartDist() 
{
    for (int i = 0; i < this->nAtoms(); i++)
    {
        for (int j = i + 1; j < this->nAtoms(); j++)
        {
            _interPartDist.push_back(coord3d::dist(_coordinates[i], _coordinates[j]));
        }
    }
    sort(_interPartDist.begin(), _interPartDist.end());
}


void structure::propertyAdjMatrix (vector<double> &p)
{
    vector< vector<int> > adjMatrix = structure::createAdjMatrix(p);
    _adjMatrix = adjMatrix;
}

void structure::propertyDistMatrix()
{
    for (int i = 0; i < this->nAtoms(); i++)
    {
        vector<double> currentRow;
        for (int j = 0; j < this->nAtoms(); j++)
        {
            currentRow.push_back(coord3d::dist(_coordinates[i], _coordinates[j]));
        }
        _distMatrix.push_back(currentRow);
    }
}

void structure::propertyGraph(double rm, double eps)
{
    undirectedGraph G = this->createGraph(rm, eps);
    _uGraph = G;
}


void structure::propertyGraph_ignoreCenter(double rm, double eps)
{
    undirectedGraph G = this->createGraph_ignoreCenter(rm, eps);
    _uGraph = G;
}

double structure::longest_nearest_neighbour_distance (unsigned int n)
{
    double distance(0);

    
    for (vector<coord3d>::size_type i = 0; i < _coordinates.size(); i++)
    {
        vector<double> allDistances;
        for (vector<coord3d>::size_type j = 0; j < _coordinates.size(); j++)
        {
            if (i == j) continue;
            allDistances.push_back(coord3d::dist(_coordinates[i],_coordinates[j]));
        }
        sort (allDistances.begin(), allDistances.end());
        allDistances.resize(n);
        //for (vector<double>::size_type a = 0; a < allDistances.size(); a++)
        //{
        //    cout << allDistances[a] << endl;
        //}

        if (allDistances.back() > distance)
        {
            distance = allDistances.back();
        }
    }


    return distance;
}


//symmetry
//sig mirror; i: 0=yz 1=xz, 2=xy
vector<coord3d> structure::sig (int a) {
    vector<coord3d> coord = _coordinates;
    for (int i = 0; i < this->nAtoms(); i++) { 
        coord[i][a] = - _coordinates[i][a];
    }
    return coord;
}

//c2 rotation; ab: 12=x, 02=y, 01=z
vector<coord3d> structure::c2 (int a, int b) {
    vector<coord3d> coord = _coordinates;
    for (int i = 0; i < this->nAtoms(); i++) { 
        coord[i][a] = -_coordinates[i][a];
        coord[i][b] = -_coordinates[i][b];
    }
    return coord;
}

//inversion
vector<coord3d> structure::inv () {
    vector<coord3d> coord = _coordinates;
    for (int i = 0; i < this->nAtoms(); i++) { 
        coord[i] = - _coordinates[i];
    }
    return coord;
}

void structure::randomize(int number)
{
        vector<coord3d> coordinates;
        random_device r;
        default_random_engine generator{r()};
        for (int i = 0; i < number; i++)
        {
            double x[3];
            for (int j = 0; j < 3; j++)
            {
                uniform_real_distribution<double> distribution(-1.0,1.0);
                x[j] = distribution(generator);
            }
            coord3d coord(x);
            coordinates.push_back(coord);
        }
        this->setCoordinates(coordinates);
}
