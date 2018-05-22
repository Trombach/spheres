#ifndef GRAPH
#define GRAPH

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include "geometry.h"

struct vertex_info
{
    unsigned int index;
    coord3d coordinate;
};

typedef std::pair<double,double> GraphCoords;

typedef boost::adjacency_list<  boost::vecS, 
                                boost::vecS, 
                                boost::undirectedS, 
                                //boost::property<boost::vertex_index_t, int>,
                                //boost::property<boost::vertex_coordinate_t, coord3d>,
                                vertex_info,
                                boost::property<boost::edge_index_t, int> > undirectedGraph;

template <class T> struct output_visitor : public boost::planar_face_traversal_visitor 
{
    unsigned int _F;
    unsigned int _F3, _F4, _F5, _F6, _F7, _F8, _F9, _F10, _F11, _F12;
    T &_collector;

    output_visitor(T &collector) : _F(0), _F3(0), _F4(0), _F5(0), _F6(0), _F7(0), _F8(0), _F9(0), _F10(0), _F11(0), _F12(0), _collector(collector) {}

    void begin_face() { _F = 0; }

    template <typename Vertex> void next_vertex(Vertex v)
    { 
        _F++;
    }

    void end_face() 
    {
        switch (_F)
        {
            case 3: _F3++;
                    break;
            case 4: _F4++;
                    break;
            case 5: _F5++;
                    break;
            case 6: _F6++;
                    break;
            case 7: _F7++;
                    break;
            case 8: _F8++;
                    break;
            case 9: _F9++;
                    break;
            case 10: _F10++;
                     break;
            case 11: _F11++;
                     break;
            case 12: _F12++;
                     break;
            default: std::cerr << "Face of degree " << _F << " ignored." << std::endl;
                     break;
        }
        _F = 0;
    }

    void end_traversal()
    {
        _collector.setFaces(_F3,_F4,_F5,_F6,_F7,_F8,_F9,_F10,_F11,_F12);
        //std::cout << _F3 << " " << _F4 << " " << _F5 << " " << _F6 << " " << _F7 << " " << _F8 << " " << _F9 << " " << _F10 << " " << _F11 << " " << _F12 << std::endl;
    }
};

template <typename GraphFirst, typename GraphSecond> struct print_callback 
{
    print_callback( const GraphFirst &graph1, 
                    const GraphSecond &graph2, 
                    std::vector<int> &allMatches, 
                    std::vector<std::vector<int> > &mapping) :  m_graph1(graph1),
                                                                m_graph2(graph2),
                                                                _allMatches(allMatches),
                                                                _mapping(mapping)
    {}
    template <typename CorrespondenceMapFirstToSecond, typename CorrespondenceMapSecondToFirst>
    bool operator() (   CorrespondenceMapFirstToSecond correspondence_map_1_to_2,
                        CorrespondenceMapSecondToFirst correspondence_map_2_to_1)
                        //typename boost::graph_traits<GraphFirst>::vertices_size_type subgraph_size)
    {
    unsigned int match(0);
    std::vector<int> single_mapping(12, 0);
    BGL_FORALL_VERTICES_T(vertex1, m_graph1, GraphFirst)
    {
        if (get(correspondence_map_1_to_2, vertex1) != boost::graph_traits<GraphSecond>::null_vertex())
        {
            //cout << vertex1 << " <-> " << get(correspondence_map_1_to_2, vertex1) << endl;
            single_mapping[vertex1] = get(correspondence_map_1_to_2, vertex1);
            match++;
        }
    }
    _mapping.push_back(single_mapping);
    _allMatches.push_back(match);
    //cout << "---" << endl;

    return (true);
}
                        
    private:
        std::vector<int> &_allMatches;
        std::vector<std::vector<int> > &_mapping;
        const GraphFirst &m_graph1;
        const GraphSecond &m_graph2;
};

class pos_writer 
{
    public:
        pos_writer(std::vector<GraphCoords> pos) : _positions(pos)
        {}
        template <class Vertex> 
            void operator()(std::ostream &out, const Vertex &v) const
            {
                out << " [shape=circle, pos=\"" << _positions[v].first << "," << _positions[v].second << "!\"]";
            }
    private:
        std::vector<GraphCoords> _positions;
};

#endif
