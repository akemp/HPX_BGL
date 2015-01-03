#ifndef HEADERS_BGL
#define HEADERS_BGL
#ifndef NO_HPX
#include <hpx/hpx_main.hpp>
#include <hpx/include/runtime.hpp>
#include <hpx/include/thread_executors.hpp>
#include <hpx/util/high_resolution_timer.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/iostreams.hpp>
#endif


#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <mutex>
#include <thread>         // std::thread
#include "../metis/include/metis.h"

#include <boost/random.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>


typedef std::pair < int, int > Edge;
typedef std::vector<Edge> Edges;

struct multi_name_t {
	typedef boost::vertex_property_tag kind;
};

typedef boost::property<multi_name_t, std::vector<int> > MultiColor; //parent, color, partition, distance
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	MultiColor> MultiGraph;

int get_parts(std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, std::vector<idx_t> &part,
	idx_t nparts);
void toCSR(const std::vector<std::vector<idx_t>>& nodes, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy);
#endif