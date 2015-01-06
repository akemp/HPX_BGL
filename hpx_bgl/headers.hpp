#ifndef HEADERS_BGL
#define HEADERS_BGL

#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <mutex>
#include <thread>         // std::thread
#include <future>
#include "../metis/include/metis.h"

#include <boost/random.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include "high_resolution_timer.hpp"


typedef std::pair < int, int > Edge;
typedef std::vector<Edge> Edges;

struct multi_name_t {
	typedef boost::vertex_property_tag kind;
};

typedef boost::property<multi_name_t, std::vector<int> > MultiColor; //parent, color, partition, distance
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	MultiColor> Graph;

#endif