#include <iostream>
#include <vector>
#include <iomanip>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include "../iop.h"
#include "../geometry.h"
#include "../structure.h"
#include "../graph.h"

using namespace std;
using namespace boost;


struct DegreeCollection
{
    unsigned int v2, v3, v4, v5, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;
    unsigned int number, Nc;
    undirectedGraph printGraph;
    undirectedGraph removedEdgesGraph;
    vector<int> mapping;

    bool isVertexEqual(DegreeCollection D) const
    {
        if (this->v2 == D.v2 && 
            this->v3 == D.v3 && 
            this->v4 == D.v4 && 
            this->v5 == D.v5) return true;
        return false;
    }

    bool isFaceEqual(DegreeCollection D) const
    {
        if (this->f3 == D.f3 &&
            this->f4 == D.f4 &&
            this->f5 == D.f5 &&
            this->f6 == D.f6 &&
            this->f7 == D.f7 &&
            this->f8 == D.f8 &&
            this->f9 == D.f9 &&
            this->f10 == D.f10 &&
            this->f11 == D.f11 &&
            this->f12 == D.f12) return true;
        return false;
    }

    void setFaces(  unsigned int F3,
                    unsigned int F4,
                    unsigned int F5,
                    unsigned int F6,
                    unsigned int F7,
                    unsigned int F8,
                    unsigned int F9,
                    unsigned int F10,
                    unsigned int F11,
                    unsigned int F12)
    {
        this->f3 = F3;
        this->f4 = F4;
        this->f5 = F5;
        this->f6 = F6;
        this->f7 = F7;
        this->f8 = F8;
        this->f9 = F9;
        this->f10 = F10;
        this->f11 = F11;
        this->f12 = F12;
    }

    void setVertices(   unsigned int V2,
                        unsigned int V3,
                        unsigned int V4,
                        unsigned int V5)
    {
        this->v2 = V2;
        this->v3 = V3;
        this->v4 = V4;
        this->v5 = V5;
        this->Nc = (2 * v2 + 3 * v3 + 4 * v4 + 5 * v5) / 2;
    }

};



int main (int argc, char *argv[])
{

    cout << endl;

    string fileName = argv[argc - 1];

    if (fexists(fileName))
    {
        cout << "\tFile " << fileName << " exists." << endl;
    }
    else
    {
        cout << "\tFile " << fileName << " does not exist." << endl;
        return 1;
    }

    cout << endl;

    cout << "\t+ Reading and processing structures" << endl;
    vector<structure> KS = readallstruct(fileName);
    cout << "\t\t#structures :" << KS.size() << endl;

    for (vector<structure>::size_type i = 0; i < KS.size(); i++)
    {
        KS[i].propertyGraph_ignoreCenter();
    }


    //structure of icosahedron
    const double phi = (1 + sqrt(5)) / 2;
    vector<coord3d> ico_coord =
    {
        coord3d(0, 1, phi),
        coord3d(0, 1, -phi),
        coord3d(0, -1, phi),
        coord3d(0, -1, -phi),
        coord3d(1, phi, 0),
        coord3d(1, -phi, 0),
        coord3d(-1, phi, 0),
        coord3d(-1, -phi, 0),
        coord3d(phi, 0, 1),
        coord3d(phi, 0, -1),
        coord3d(-phi, 0, 1),
        coord3d(-phi, 0, -1)
    };

    structure ico (0, ico_coord);
    ico.propertyGraph(2);


    
    //comparison
    auto isSubgraphIco = [&] (vector<int> M)
    {
        for (vector<int>::size_type i; i < M.size(); i++)
        {
            if (M[i] != 12) return false;
        }
        return true;
    };

    undirectedGraph graph2 = ico.getGraph();
    ofstream icoOut;
    icoOut.open("output/graphs/graph_ico_test.dot");
    write_graphviz(icoOut, graph2);//, pos_writer(graph2_coords));
    icoOut.close();
    

    vector<DegreeCollection> allGraphs;
    for (vector<structure>::size_type i = 0; i < KS.size(); i++)
    {
        cout << "(" << KS[i].getNumber() << ")" << " ";
        undirectedGraph graph1 = KS[i].getGraph();

        vector<int> matches;
        vector<int> mapping(12, 0);//mapping from vector index to value at that index
        print_callback<undirectedGraph, undirectedGraph> my_callback(graph1, graph2, matches, mapping);
        vf2_subgraph_mono(graph1, graph2, my_callback);
        if (isSubgraphIco(matches))
        {
            cout << "yes, " << matches.size() << " mappings" << endl;
        }
        else cout << "no" << endl;


        //create graph matching the icosahedral graph
        undirectedGraph print_graph1(12);
        BGL_FORALL_EDGES_T(edge, graph1, undirectedGraph)
        {
            //cout << source(edge, graph1) << " - " << target(edge, graph1) << endl;
            add_edge(mapping[source(edge, graph1)], mapping[target(edge, graph1)], print_graph1);
        }

        //add coordinate properties to print_graph1
        BGL_FORALL_VERTICES_T(vertex, graph1, undirectedGraph)
        {
            unsigned int index = graph1[vertex].index;
            coord3d coord = graph1[vertex].coordinate;

            print_graph1[mapping[index]].coordinate = coord;
        }

        //create graph of removed edges collection by removing edges from icosahedral graph
        undirectedGraph removed_edges_graph = graph2;
        BGL_FORALL_EDGES_T(edge, print_graph1, undirectedGraph)
        {
            remove_edge(source(edge, print_graph1), target(edge, print_graph1), removed_edges_graph);
        }

        //save mapped graph to DegreeCollection for later use
        
        DegreeCollection graphDegrees;
        graphDegrees.number = KS[i].getNumber();
        graphDegrees.printGraph = print_graph1;
        graphDegrees.removedEdgesGraph = removed_edges_graph;
        graphDegrees.mapping = mapping;

        ofstream graphOut;
        graphOut.open("output/graphs/graph" 
                + to_string(KS[i].getNumber()) + "-" + to_string(num_edges(print_graph1)) + ".dot");
        write_graphviz(graphOut, print_graph1);
        graphOut.close();

        //counting faces and vertices

        output_visitor<DegreeCollection> my_visitor(graphDegrees);

        //initialize edge index
        property_map<undirectedGraph, edge_index_t>::type e_index = get(edge_index, graph1);
        graph_traits<undirectedGraph>::edges_size_type edge_count = 0;
        graph_traits<undirectedGraph>::edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(graph1); ei != ei_end; ++ei)
        {
            put(e_index, *ei, edge_count++);
        }

        typedef vector<graph_traits<undirectedGraph>::edge_descriptor > vec_t;
        vector<vec_t> embedding(num_vertices(graph1));
        boyer_myrvold_planarity_test(   boyer_myrvold_params::graph = graph1,
                                        boyer_myrvold_params::embedding = &embedding[0]);

        planar_face_traversal(graph1, &embedding[0], my_visitor);

        unsigned int v2(0), v3(0), v4(0), v5(0);
        BGL_FORALL_VERTICES_T(vertex, graph1, undirectedGraph)
        {
            if (degree(vertex, graph1) == 2) v2++;
            if (degree(vertex, graph1) == 3) v3++;
            if (degree(vertex, graph1) == 4) v4++;
            if (degree(vertex, graph1) == 5) v5++;
        }
        graphDegrees.setVertices(v2,v3,v4,v5);
        
        allGraphs.push_back(graphDegrees);
    }

    //analysis

    vector< vector< vector< DegreeCollection > > > vertex_face_Degrees_equality_classes;
    //sort by vertices first
    for (vector<DegreeCollection>::size_type i = 0; i < allGraphs.size(); i++)
    {
        if (vertex_face_Degrees_equality_classes.size() == 0)
        {
            vector<DegreeCollection> face_eqClass;
            vector< vector< DegreeCollection > > vertex_eqClass;
            face_eqClass.push_back(allGraphs[i]);
            vertex_eqClass.push_back(face_eqClass);
            vertex_face_Degrees_equality_classes.push_back(vertex_eqClass);
        }
        else
        {
            bool matched_vertices(false);
            bool matched_faces(false);
            for (   vector< vector< vector< DegreeCollection > > >::size_type j = 0; 
                    j < vertex_face_Degrees_equality_classes.size(); 
                    j++)
            {
                for (   vector< vector< DegreeCollection > >::size_type k = 0; 
                        k < vertex_face_Degrees_equality_classes[j].size(); 
                        k++)
                {
                    if (allGraphs[i].isVertexEqual(vertex_face_Degrees_equality_classes[j][k][0]))
                    {
                        matched_vertices = true;
                        if (allGraphs[i].isFaceEqual(vertex_face_Degrees_equality_classes[j][k][0]))
                        {
                            matched_faces = true;
                            vertex_face_Degrees_equality_classes[j][k].push_back(allGraphs[i]);
                            //cout << i << " " << j << " " << k << " 1" << endl;
                            break;
                        }
                    }
                }
                if (matched_vertices && matched_faces == false)
                {
                    vector<DegreeCollection> face_newEqClass;
                    face_newEqClass.push_back(allGraphs[i]);
                    vertex_face_Degrees_equality_classes[j].push_back(face_newEqClass);
                    //cout << i << " " << j << " 2" << endl;
                    break;
                }
            }
            if (matched_vertices == false && matched_faces == false)
            {
                vector<DegreeCollection> A;
                vector< vector< DegreeCollection > > B;
                A.push_back(allGraphs[i]);
                B.push_back(A);
                vertex_face_Degrees_equality_classes.push_back(B);
                //cout << i << " 3" << endl;
            }
            //cout << i << " " << matched_vertices << " " << matched_faces << endl;
        }
    }

    auto compare_vertex = [&] (vector< vector< DegreeCollection > > a, vector< vector< DegreeCollection > > b)
    {
       if (a[0][0].Nc > b[0][0].Nc) return true; 
       if (a[0][0].Nc == b[0][0].Nc && a[0][0].v5 > b[0][0].v5) return true; 
       if (a[0][0].Nc == b[0][0].Nc && a[0][0].v5 == b[0][0].v5 && a[0][0].v4 > b[0][0].v4) return true; 
       if (a[0][0].Nc == b[0][0].Nc && a[0][0].v5 == b[0][0].v5 && a[0][0].v4 == b[0][0].v4 && a[0][0].v3 > b[0][0].v3) return true; 
       if (a[0][0].Nc == b[0][0].Nc && a[0][0].v5 == b[0][0].v5 && a[0][0].v4 == b[0][0].v4 && a[0][0].v3 == b[0][0].v3 && a[0][0].v2 > a[0][0].v2) return true; 
       return false;
    };

    auto compare_faces = [&] (DegreeCollection a, DegreeCollection b)
    {
        if (a.f3 > b.f3) return true;
        if (a.f3 == b.f3 && a.f4 > b.f4) return true;
        if (a.f3 == b.f3 && a.f4 == b.f4 && a.f5 > b.f5) return true; 
        if (a.f3 == b.f3 && a.f4 == b.f4 && a.f5 == b.f5 && a.f6 > b.f6) return true; 
        if (a.f3 == b.f3 && a.f4 == b.f4 && a.f5 == b.f5 && a.f6 == b.f6 && a.f7 > b.f7) return true; 
        if (a.f3 == b.f3 && a.f4 == b.f4 && a.f5 == b.f5 && a.f6 == b.f6 && a.f7 == b.f7 && a.f8 > b.f8) return true; 
        if (a.f3 == b.f3 && a.f4 == b.f4 && a.f5 == b.f5 && a.f6 == b.f6 && a.f7 == b.f7 && a.f8 == b.f8 && a.f9 > b.f9) return true; 
        if (a.f3 == b.f3 && a.f4 == b.f4 && a.f5 == b.f5 && a.f6 == b.f6 && a.f7 == b.f7 && a.f8 == b.f8 && a.f9 == b.f9 && a.f10 > b.f10) return true; 
        if (a.f3 == b.f3 && a.f4 == b.f4 && a.f5 == b.f5 && a.f6 == b.f6 && a.f7 == b.f7 && a.f8 == b.f8 && a.f9 == b.f9 && a.f10 == b.f10 && a.f11 > b.f11) return true; 
        if (a.f3 == b.f3 && a.f4 == b.f4 && a.f5 == b.f5 && a.f6 == b.f6 && a.f7 == b.f7 && a.f8 == b.f8 && a.f9 == b.f9 && a.f10 == b.f10 && a.f11 == b.f11 && a.f12 > b.f12) return true; 
        return false;
    };

    sort(   vertex_face_Degrees_equality_classes.begin(), 
            vertex_face_Degrees_equality_classes.end(), 
            compare_vertex);

    for (   vector< vector< vector< DegreeCollection > > >::size_type i = 0; 
            i < vertex_face_Degrees_equality_classes.size(); 
            i++)
    {
        sort(   vertex_face_Degrees_equality_classes[i][0].begin(), 
                vertex_face_Degrees_equality_classes[i][0].end(), 
                compare_faces);
    }

    for (   vector< vector< vector< DegreeCollection > > >::size_type i = 0; 
            i < vertex_face_Degrees_equality_classes.size(); 
            i++)
    {
        for (   vector< vector< DegreeCollection > >::size_type j = 0; 
                j < vertex_face_Degrees_equality_classes[i].size(); 
                j++)
        {
            cout << "Nc " << vertex_face_Degrees_equality_classes[i][j][0].Nc << " "
                << "V " << vertex_face_Degrees_equality_classes[i][j][0].v5 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].v4 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].v3 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].v2 << " F "
                << vertex_face_Degrees_equality_classes[i][j][0].f3 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].f4 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].f5 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].f6 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].f7 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].f8 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].f9 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].f10 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].f11 << " "
                << vertex_face_Degrees_equality_classes[i][j][0].f12 << " " 
                << vertex_face_Degrees_equality_classes[i][j].size() << endl;
        }
    }

    //print output files
    unsigned int num(1);
    for (   vector< vector< vector< DegreeCollection > > >::size_type i = 0; 
            i < vertex_face_Degrees_equality_classes.size(); 
            i++)
    {
        for (   vector< vector< DegreeCollection > >::size_type j = 0; 
                j < vertex_face_Degrees_equality_classes[i].size(); 
                j++)
        {
            for (   vector< DegreeCollection >::size_type k = 0; 
                    k < vertex_face_Degrees_equality_classes[i][j].size(); 
                    k++)
            {
                int Nc = vertex_face_Degrees_equality_classes[i][j][k].Nc;

                undirectedGraph print_graph1 = 
                    vertex_face_Degrees_equality_classes[i][j][k].printGraph;
    
                //output graph file
                ofstream sorted_graphOut;
                sorted_graphOut.open("output/sorted-graphs/graph" 
                        + to_string(num) + "-" + to_string(Nc) + ".dot");
                write_graphviz(sorted_graphOut, print_graph1);
                sorted_graphOut.close();


                //output coord file
                ofstream sorted_coordOut;
                sorted_coordOut.open("output/coords/coord" 
                        + to_string(num) + "-" + to_string(Nc) + ".xyz");
                int number = vertex_face_Degrees_equality_classes[i][j][k].number;
                xyzout(KS[number - 1], "coords/coord" 
                        + to_string(num) + "-" + to_string(Nc) + ".xyz");
                

                //output list of removed edges
                undirectedGraph removed_edges_graph 
                    = vertex_face_Degrees_equality_classes[i][j][k].removedEdgesGraph;

                cout << num << " " << vertex_face_Degrees_equality_classes[i][j][k].Nc << " ";
                BGL_FORALL_EDGES_T(edge, removed_edges_graph, undirectedGraph)
                {
                    cout << "(" << source(edge, removed_edges_graph) << "," 
                        << target(edge, removed_edges_graph) << ")";

                }

                double distance(0);
                BGL_FORALL_EDGES_T(edge, removed_edges_graph, undirectedGraph)
                {
                    double dist 
                        = coord3d::dist(print_graph1[source(edge, removed_edges_graph)].coordinate,
                                        print_graph1[target(edge, removed_edges_graph)].coordinate);

                if (dist > distance) distance = dist;
                }
                cout << " " << distance;

                //double nnDistance = KS[number - 1].longest_nearest_neighbour_distance(7);
                //cout << " " << nnDistance << endl;
                cout << endl;

                num++;
            }
        }
    }


                    



    return 0;
}
