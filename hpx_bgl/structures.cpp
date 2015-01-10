//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#define NO_std
#include "headers.hpp"

int get_parts(std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, std::vector<idx_t> &part,
	idx_t nparts)
{
	// Needed by parmetis
	idx_t nvtxs = xadj.size() - 1, ncon = 1;
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
void toCSR(const std::vector<std::vector<idx_t>>& nodes, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy)
{
	int total = 0;
	xadj.push_back(0);
	for (int i = 0; i < nodes.size(); ++i)
	{
		for (int j = 0; j < nodes[i].size(); ++j)
		{
			adjncy.push_back(nodes[i][j]);
			++total;
		}
		xadj.push_back(total);
	}
	return;
}

void parallel_edge_gen(vector<int>::iterator pedges, vector<vector<int>>* nodes, int size)
{

	boost::random::mt19937 rng;
	boost::random::uniform_int_distribution<> randnodes(0, nodes->size() - 1);
	rng.seed(*pedges);
	for (int i = 0; i < size; ++i)
	{
		int v0 = randnodes(rng);//pedges->v0_low;
		int v1 = randnodes(rng);//pedges->v1_low;
		if (v0 == v1)
			continue;
		{
			//undirected so no changes to final edgelist
			if (v1 < v0)
			{
				int temp = v0;
				v0 = v1;
				v1 = temp;
			}
			if (std::find((*nodes)[v0].begin(), (*nodes)[v0].end(), v1) != (*nodes)[v0].end())
				continue;
			(*nodes)[v0].push_back(v1);
		}
	}
}