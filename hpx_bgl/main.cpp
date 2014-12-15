//  Copyright (c) 2011 Matthew Anderson
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include "../metis/include/metis.h"
//#include <hpx/hpx.hpp>
//#include <hpx/hpx_init.hpp>

#include <hpx/util/high_resolution_timer.hpp>

#include <boost/throw_exception.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>


typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS > graph_t;

//int64_t *deg, *next;
typedef std::pair < int, int > Edge;
typedef std::vector<Edge> Edges;

int getres(std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, std::vector<idx_t> &part,
	idx_t nparts)
{

	// Needed by parmetis
	idx_t nvtxs = xadj.size()-1, ncon = 1;
	idx_t *vsize = NULL;
	idx_t *vwgt = NULL, *adjwgt = NULL;

	real_t *tpwgts = NULL, *ubvec = NULL;
	idx_t *options = NULL;
	idx_t objval;
	int result = METIS_PartGraphKway(
		&nvtxs, &ncon, &xadj[0], &adjncy[0], vwgt,
		vsize, adjwgt, &nparts, tpwgts,
		ubvec, options, &objval, &part[0]);
	return result;
}

void toCSR(const std::vector<std::vector<idx_t>>& links, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy)
{
	int total = 0;
	xadj.push_back(0);
	for (int i = 0; i < links.size(); ++i)
	{
		for (int j = 0; j < links[i].size(); ++j)
		{
			adjncy.push_back(links[i][j]);
			++total;
		}
		xadj.push_back(total);
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////
int init()//hpx_main(boost::program_options::variables_map &vm)
{
	using namespace std;
	int result;
	//from librelocus

	vector<idx_t> xadj;
	vector<idx_t> adjncy;

	vector<vector<idx_t>> edges(32);
	for (int i = 0; i < edges.size(); ++i)
	{
		vector<idx_t> edger;
		for (int j = 0; j < rand() % 6; ++j)
		{
			edger.push_back(rand() % edges.size());
		}
		edges[i] = edger;
	}
	toCSR(edges, xadj, adjncy);

	vector<idx_t> part(xadj.size());
	int res = getres(xadj, adjncy, part, 3);
	cout << res << endl;
	cout << "npart is " << '\n';
	for (int i = 0; i<part.size(); i++){
		cout << part[i] << " ";
	}
	cout << '\n';

//	delete xadj;
//	delete adjncy;

	cout << "Metis returned successfully." << '\n';

	return 0;
}
int main() {
    return init();//hpx::init();
}