
#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>
#include <boost/geometry/multi/geometries/multi_linestring.hpp>

using namespace std;
using namespace boost;


const float PI = 3.14159265359f;
const float HALF_PI = 1.57079632679f;

typedef boost::geometry::model::d2::point_xy<double> Point;
typedef boost::geometry::model::linestring<Point> Line;
typedef boost::geometry::model::multi_linestring <Line> multilinestring;
typedef boost::geometry::model::polygon<Point> polygon;
typedef boost::geometry::model::multi_polygon <polygon> multipolygon;

template <class T>
const std::vector<
    boost::geometry::model::d2::point_xy<T>
>
GetLineLineIntersections(
const boost::geometry::model::linestring<
boost::geometry::model::d2::point_xy<T>
> line1,
const boost::geometry::model::linestring<
boost::geometry::model::d2::point_xy<T>
> line2)
{

    std::vector<Point> points;
    boost::geometry::intersection(line1, line2, points);
    //assert(points.empty() || points.size() == 1);
    return points;
}

template <class T>
const boost::geometry::model::linestring<boost::geometry::model::d2::point_xy<T>
>
CreateLine(const std::vector<boost::geometry::model::d2::point_xy<T> >& v)
{
    return boost::geometry::model::linestring<
        boost::geometry::model::d2::point_xy<T>
    >(v.begin(), v.end());
}

struct fuzzy_equal_to
    : public std::binary_function<double, double, bool>
{
    fuzzy_equal_to(const double tolerance = std::numeric_limits<double>::epsilon())
        : m_tolerance(tolerance)
    {
        assert(tolerance >= 0.0);
    }
    bool operator()(const double lhs, const double rhs) const
    {
        return rhs > (1.0 - m_tolerance) * lhs
            && rhs < (1.0 + m_tolerance) * lhs;
    }
    const double m_tolerance;
};

void splitLine(vector<Line> spotsl, vector<Line>& roadsegs)
{
    for (int i = 0; i < spotsl.size(); ++i)
    {
        Line line1 = spotsl[i];
        Line seg;
        vector<vector<Point>> pts;
        for (int j = 0; j < spotsl.size(); ++j)
        {
            Line line2 = spotsl[j];
            vector<Point> ints = GetLineLineIntersections(line1, line2);
            if (ints.size() > 0)
            {
                vector<Point> temp;
                temp.push_back(ints[0]);
                temp.push_back(spotsl[i].front());
                pts.push_back(temp);
            }
        }
        sort(pts.begin(), pts.end(), [](const vector<Point>& p1, const vector<Point>& p2)
        {return boost::geometry::distance(p1.front(), p1.back()) < boost::geometry::distance(p2.front(), p2.back()); });
        for (int j = 0; j < pts.size(); ++j)
        {
            seg.push_back(pts[j][0]);
        }
        roadsegs.push_back(seg);
    }
    return;
}

typedef adjacency_list < vecS, vecS, undirectedS > graph_t;

template < typename TimeMap > class bfs_time_visitor :public default_bfs_visitor {
    typedef typename property_traits < TimeMap >::value_type T;
public:
    bfs_time_visitor(TimeMap tmap, T & t) :m_timemap(tmap), m_time(t) { }
    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph & g) const
    {
        put(m_timemap, u, m_time++);
        if (transfermap[u] != index)
        {
            
        }
    }
    int index;
    vector<int> transfermap;

    TimeMap m_timemap;
    T & m_time;
};
void search()
{
    /*
    using namespace boost;
    // Select the graph type we wish to use
    // Set up the vertex IDs and names
    // Specify the edges in the graph
    typedef std::pair < int, int >E;
    graph_t g;
    vector<E> edges;

    vector<Line> spots;
    for (int i = 0; i < edges.size(); ++i)
    {
        Line l;
        l.push_back((points[e.first]));
        l.push_back(points[e.second]));
        spots.push_back(l);
    }
    vector<Line> segs;
    splitLine(spots, segs);
    segs = spots;

    // Typedefs
    typedef graph_traits < graph_t >::vertices_size_type Size;

    // a vector to hold the discover time property for each vertex
    std::vector < Size > dtime(num_vertices(g));
    typedef
        iterator_property_map<std::vector<Size>::iterator,
        property_map<graph_t, vertex_index_t>::const_type>
        dtime_pm_type;
    dtime_pm_type dtime_pm(dtime.begin(), get(vertex_index, g));

    Size time = 0;
    bfs_time_visitor < dtime_pm_type >vis(dtime_pm, time);
    breadth_first_search(g, vertex(0, g), visitor(vis));

    */
}