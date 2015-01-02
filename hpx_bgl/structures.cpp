//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define NO_HPX
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