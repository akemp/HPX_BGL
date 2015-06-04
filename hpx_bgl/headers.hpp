#ifndef HEADERS_BGL
#define HEADERS_BGL

#include <hpx/hpx_main.hpp>
#include <hpx/util/high_resolution_timer.hpp>
#include <hpx/include/components.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include "../metis/include/metis.h"

#include <boost/random.hpp>
#include <boost/graph/adjacency_list.hpp>


typedef std::pair < int, int > Edge;
typedef std::vector<Edge> Edges;


void parallel_edge_gen(std::vector<int>::iterator pedges, std::vector<std::vector<int>>& nodes, int size);
void toCSR(const std::vector<std::vector<idx_t>>& nodes, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy);
int get_parts(std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, std::vector<idx_t> &part,
	idx_t nparts);
#endif