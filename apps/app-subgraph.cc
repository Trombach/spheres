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

    //coordinate pairs for ico graph
    vector<GraphCoords> graph2_coords =
    {
        make_pair(0,0),
        make_pair(500,454.6),
        make_pair(500,115.5),
        make_pair(573.5,331.435),
        make_pair(500,866),
        make_pair(650,201.435),
        make_pair(350,370),
        make_pair(500,240),
        make_pair(1000,0),
        make_pair(650,370),
        make_pair(350,201.435),
        make_pair(426.5,331.435)
    };

    
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
    
    for (vector<structure>::size_type i = 0; i < KS.size(); i++)
    {
        cout << "(" << KS[i].getNumber() << ")" << " ";
        undirectedGraph graph1 = KS[i].getGraph();

        vector<int> matches;
        vector<int> mapping(12, 0);
        print_callback<undirectedGraph, undirectedGraph> my_callback(graph1, graph2, matches, mapping);
        vf2_subgraph_mono(graph1, graph2, my_callback);
        if (isSubgraphIco(matches))
        {
            cout << "yes, " << matches.size() << " mappings" << endl;
        }
        else cout << "no" << endl;

        vector<GraphCoords> graph1_coords(12);
        for (vector<int>::size_type j = 0; j < mapping.size(); j++)
        {
            //cout << j << " - " << mapping[j] << endl;
            graph1_coords[mapping[j]] = graph2_coords[j];
        }

        undirectedGraph print_graph1(12);
        BGL_FORALL_EDGES_T(edge, graph1, undirectedGraph)
        {
            //cout << source(edge, graph1) << " - " << target(edge, graph1) << endl;
            add_edge(mapping[source(edge, graph1)], mapping[target(edge, graph1)], print_graph1);
        }


        ofstream graphOut;
        graphOut.open("output/graphs/graph" + to_string(KS[i].getNumber()) + "-" + to_string(num_edges(print_graph1)) + ".dot");
        write_graphviz(graphOut, print_graph1);//, pos_writer(graph1_coords));
        graphOut.close();

        output_visitor my_visitor;

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

    }


    return 0;
}
