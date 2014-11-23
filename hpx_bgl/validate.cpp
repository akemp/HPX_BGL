//  Copyright (c) 2011 Matthew Anderson
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <boost/assert.hpp>

#include "../graph500/make-edgelist.h"
#include "../graph500/graph500.h"
#include "../graph500/generator/make_graph.h"

#include "../graph500/prng.h"
#include "../graph500/generator/splittable_mrg.h"

#define VERBOSE true

static int compare_doubles(const void* a, const void* b) {
  double aa = *(const double*)a;
  double bb = *(const double*)b;
  return (aa < bb) ? -1 : (aa == bb) ? 0 : 1;
}

void get_maxedge(struct packed_edge *IJ_in, int64_t nedge, int64_t &maxvtx)
{
	packed_edge * IJ = IJ_in;
	//nIJ = nedge;
	maxvtx = -1;

	int64_t k;
	for (k = 0; k < nedge; ++k) {
		if (get_v0_from_edge(&IJ[k]) > maxvtx)
			maxvtx = get_v0_from_edge(&IJ[k]);
		if (get_v1_from_edge(&IJ[k]) > maxvtx)
			maxvtx = get_v1_from_edge(&IJ[k]);
	}
}

int create_graph_from_edgelist(struct packed_edge *IJ_in, int64_t nedge,
	int64_t &maxvtx, int64_t& maxdeg, int64_t* head, int64_t* deg, int64_t* next)
{
	int err = 0;

	maxdeg = -1;
	packed_edge * IJ = IJ_in;
	int64_t k;
	if (!head) return -1;
	for (k = 0; k < nedge; ++k) {
		const int64_t i = get_v0_from_edge(&IJ[k]);
		const int64_t j = get_v1_from_edge(&IJ[k]);
		int64_t t_head, t;

		if (i >= 0 && j >= 0 && i != j) {
			next[2 * k] = -1;
			next[1 + 2 * k] = -1;
			t = 2 * k + 1; /* Point at the *other* end. */
			t_head = head[i];
			head[i] = t;
			assert(t_head < 2 * nedge);
			next[t] = t_head;
			++deg[i];

			--t;
			t_head = head[j];
			head[j] = t;
			assert(t_head < 2 * nedge);
			next[t] = t_head;
			++deg[j];
		}
	}

	for (int64_t kg = 0; kg <= maxvtx; ++kg)
		if (deg[kg] > maxdeg)
			maxdeg = deg[kg];

	return err;
}

int make_bfs_tree(int64_t *bfs_tree_out, int64_t *max_vtx_out,
	int64_t srcvtx, packed_edge * IJ, int64_t maxvtx, int64_t* head, int64_t* next)
{
	int64_t * bfs_tree = bfs_tree_out;
	int err = 0;
	const int64_t nv = maxvtx + 1;

	int64_t k, k1, k2, newk2;
	int64_t * vlist;

	*max_vtx_out = maxvtx;

	bfs_tree[srcvtx] = srcvtx;
	newk2 = 1;

	vlist = new int64_t[nv];
	if (!vlist) return -1;
	bfs_tree[srcvtx] = srcvtx;
	k1 = 0; k2 = 1;
	vlist[0] = srcvtx;

	for (k = 0; k < srcvtx; ++k)
		bfs_tree[k] = -1;
	for (k = srcvtx + 1; k < nv; ++k)
		bfs_tree[k] = -1;

	while (k1 != k2) {
		int64_t k, newk2 = k2;
		for (k = k1; k < k2; ++k) {
			const int64_t parent = vlist[k];
			int64_t p = head[parent];
			while (p >= 0) {
				const int64_t newv = ((p % 2) ? get_v1_from_edge(&IJ[p / 2]) : get_v0_from_edge(&IJ[p / 2]));
				if (bfs_tree[newv] < 0) {
					bfs_tree[newv] = parent;
					vlist[newk2++] = newv;
				}
				p = next[p];
			}
			k1 = k2;
			k2 = newk2;
		}
	}
	free(vlist);

	return err;
}


void verify(int NBFS_max, int64_t nvtx_scale, int NBFS, int64_t nedge,
	packed_edge* IJ, int64_t* bfs_root)
{

	int * has_adj;
	int m, err;
	int64_t k, t;

	has_adj = new int[nvtx_scale];
	for (k = 0; k < nvtx_scale; ++k)
		has_adj[k] = 0;
	for (k = 0; k < nedge; ++k) {
		const int64_t i = get_v0_from_edge(&IJ[k]);
		const int64_t j = get_v1_from_edge(&IJ[k]);
		if (i != j)
			has_adj[i] = has_adj[j] = 1;
	}

	/* Sample from {0, ..., nvtx_scale-1} without replacement. */
	m = 0;
	t = 0;
	while (m < NBFS && t < nvtx_scale) {
		//std::cout << m << std::endl;
		double R = mrg_get_double_orig(prng_state);
		if (!has_adj[t] || (nvtx_scale - t)*R > NBFS - m) ++t;
		else bfs_root[m++] = t++;
	}
	if (t >= nvtx_scale && m < NBFS) {
		if (m > 0) {
			fprintf(stderr, "Cannot find %d sample roots of non-self degree > 0, using %d.\n",
				NBFS, m);
			NBFS = m;
		}
		else {
			fprintf(stderr, "Cannot find any sample roots of non-self degree > 0.\n");
			exit(EXIT_FAILURE);
		}
	}

	free(has_adj);
}


void search(int NBFS, int64_t nvtx_scale, int64_t* bfs_root, int64_t nedge, packed_edge* IJ, int64_t maxvtx, int64_t* head, int64_t *next)
{
	//bool VERBOSE = true;
	int m, err;
	int64_t k, t;
	int64_t* bfs_nedge = new int64_t[nvtx_scale];
	for (m = 0; m < NBFS; ++m) {
		int64_t *bfs_tree, max_bfsvtx;

		/* Re-allocate. Some systems may randomize the addres... */
		bfs_tree = new int64_t[nvtx_scale];
		assert(bfs_root[m] < nvtx_scale);

		if (VERBOSE) fprintf(stderr, "Running bfs %d...", m);

		err = make_bfs_tree(bfs_tree, &max_bfsvtx, bfs_root[m], &(*IJ), maxvtx, head, next);

		if (VERBOSE) fprintf(stderr, "done\n");

		if (err) {
			perror("make_bfs_tree failed");
			abort();
		}

		if (VERBOSE) fprintf(stderr, "Verifying bfs %d...", m);
		bfs_nedge[m] = verify_bfs_tree(bfs_tree, max_bfsvtx, bfs_root[m], IJ, nedge);
		if (VERBOSE) fprintf(stderr, "done\n");
		if (bfs_nedge[m] < 0) {
			fprintf(stderr, "bfs %d from %" PRId64 " failed verification (%" PRId64 ")\n",
				m, bfs_root[m], bfs_nedge[m]);
			abort();
		}

		free(bfs_tree);
	}
}