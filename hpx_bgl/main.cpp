//  Copyright (c) 2011 Matthew Anderson
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
//#include <hpx/hpx.hpp>
//#include <hpx/hpx_init.hpp>

#include <hpx/util/high_resolution_timer.hpp>

#include "../graph500/make-edgelist.h"
#include "../graph500/graph500.h"
#include "../graph500/generator/make_graph.h"

#include "../graph500/prng.h"
#include "../graph500/generator/splittable_mrg.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>


typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS > graph_t;

//int64_t *deg, *next;
typedef std::pair < int, int > Edge;
typedef std::vector<Edge> Edges;

#define VERBOSE true

void get_maxedge(struct packed_edge *IJ_in, int64_t nedge, int64_t &maxvtx);

int create_graph_from_edgelist(struct packed_edge *IJ_in, int64_t nedge,
	int64_t &maxvtx, int64_t& maxdeg, int64_t* head, int64_t* deg, int64_t* next);

int make_bfs_tree(int64_t *bfs_tree_out, int64_t *max_vtx_out,
	int64_t srcvtx, packed_edge * IJ, int64_t maxvtx, int64_t* head, int64_t* next);

void verify(int NBFS_max, int64_t nvtx_scale, int NBFS, int64_t nedge,
	packed_edge* IJ, int64_t* bfs_root);

void search(int NBFS, int64_t nvtx_scale, int64_t* bfs_root, int64_t nedge,
	packed_edge* IJ, int64_t maxvtx, int64_t* head, int64_t *next);

void verify_boost(int NBFS_max, int64_t nvtx_scale, int NBFS, int64_t nedge,
	Edges* edges, int64_t* bfs_root)
{
	using namespace boost;
	int * has_adj;
	int m, err;
	int64_t k, t;

	has_adj = new int[nvtx_scale];
	for (k = 0; k < nvtx_scale; ++k)
		has_adj[k] = 0;

	graph_traits<graph_t>::edge_iterator ei, ei_end;
	for (auto it = edges->begin(); it < edges->end(); ++it)//(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
	{

		const int64_t i = it->first;//source(*ei, g);
		const int64_t j = it->second;//target(*ei, g);
		if (i != j)
			has_adj[i] = has_adj[j] = 1;
	}

	/* Sample from {0, ..., nvtx_scale-1} without replacement. */
	m = 0;
	t = 0;
	while (m < NBFS && t < nvtx_scale) {
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


static int
compute_levels_boost(int64_t * level,
int64_t nv, const int64_t * bfs_tree, int64_t root)
{
	int err = 0;

	{
		int terr;
		int64_t k;

			for (k = 0; k < nv; ++k)
				level[k] = (k == root ? 0 : -1);

			for (k = 0; k < nv; ++k) {
			if (level[k] >= 0) continue;
			terr = err;
			if (!terr && bfs_tree[k] >= 0 && k != root) {
				int64_t parent = k;
				int64_t nhop = 0;
				/* Run up the tree until we encounter an already-leveled vertex. */
				while (parent >= 0 && level[parent] < 0 && nhop < nv) {
					assert(parent != bfs_tree[parent]);
					parent = bfs_tree[parent];
					++nhop;
				}
				if (nhop >= nv) terr = -1; /* Cycle. */
				if (parent < 0) terr = -2; /* Ran off the end. */

				if (!terr) {
					/* Now assign levels until we meet an already-leveled vertex */
					/* NOTE: This permits benign races if parallelized. */
					nhop += level[parent];
					parent = k;
					while (level[parent] < 0) {
						assert(nhop > 0);
						level[parent] = nhop--;
						parent = bfs_tree[parent];
					}
					assert(nhop == level[parent]);

					/* Internal check to catch mistakes in races... */
#if !defined(NDEBUG)
					nhop = 0;
					parent = k;
					int64_t lastlvl = level[k] + 1;
					while (level[parent] > 0) {
						assert(lastlvl == 1 + level[parent]);
						lastlvl = level[parent];
						parent = bfs_tree[parent];
						++nhop;
					}
#endif
				}
			}
			if (terr) { err = terr;}
			}
	}
	return err;
}

int64_t
verify_bfs_tree_boost(int64_t *bfs_tree_in, int64_t max_bfsvtx,
int64_t root,
Edges* edges, int64_t nedge)
{
	int64_t * bfs_tree = bfs_tree_in;
	//const packed_edge * IJ = IJ_in;

	int err, nedge_traversed;
	int64_t * seen_edge, * level;

	const int64_t nv = max_bfsvtx + 1;

	/*
	This code is horrifically contorted because many compilers
	complain about continue, return, etc. in parallel sections.
	*/

	if (root > max_bfsvtx || bfs_tree[root] != root)
		return -999;

	err = 0;
	nedge_traversed = 0;
	seen_edge = new int64_t[2 * (nv)];
	level = &seen_edge[nv];

	err = compute_levels_boost(level, nv, bfs_tree, root);

	if (err) goto done;

	 {
		int64_t k;
		int terr = 0;
			for (k = 0; k < nv; ++k)
				seen_edge[k] = 0;

			//for (k = 0; k < nedge; k++) {
			boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
			//for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
			for (auto it = edges->begin(); it < edges->end(); ++it) {
				const int64_t i = it->first;//source(*ei, g);
				const int64_t j = it->second;//target(*ei, g);
				int64_t lvldiff;
				terr = err;

				if (i < 0 || j < 0) continue;
				if (i > max_bfsvtx && j <= max_bfsvtx) terr = -10;
				if (j > max_bfsvtx && i <= max_bfsvtx) terr = -11;
				if (terr) { err = terr; }
				if (terr || i > max_bfsvtx /* both i & j are on the same side of max_bfsvtx */)
					continue;

				/* All neighbors must be in the tree. */
				if (bfs_tree[i] >= 0 && bfs_tree[j] < 0) terr = -12;
				if (bfs_tree[j] >= 0 && bfs_tree[i] < 0) terr = -13;
				if (terr) { err = terr; }
				if (terr || bfs_tree[i] < 0 /* both i & j have the same sign */)
					continue;

				/* Both i and j are in the tree, count as a traversed edge.

				NOTE: This counts self-edges and repeated edges.  They're
				part of the input data.
				*/
				++nedge_traversed;
				/* Mark seen tree edges. */
				if (i != j) {
					if (bfs_tree[i] == j)
						seen_edge[i] = 1;
					if (bfs_tree[j] == i)
						seen_edge[j] = 1;
				}
				lvldiff = level[i] - level[j];
				/* Check that the levels differ by no more than one. */
				if (lvldiff > 1 || lvldiff < -1)
					terr = -14;
				if (terr) { err = terr; }
			}

		if (!terr) {
			/* Check that every BFS edge was seen and that there's only one root. */
				for (k = 0; k < nv; ++k) {
				terr = err;
				if (!terr && k != root) {
					if (bfs_tree[k] >= 0 && !seen_edge[k])
						terr = -15;
					if (bfs_tree[k] == k)
						terr = -16;
					if (terr) { err = terr; }
				}
				}
		}
	}
done:

	free(seen_edge);
	if (err) return err;
	return nedge_traversed;
}

int make_bfs_tree_boost(int64_t *bfs_tree_out,
	int64_t srcvtx, Edges *edges, int64_t maxvtx, int64_t* head, int64_t* next)
{
	int64_t * bfs_tree = bfs_tree_out;
	int err = 0;
	const int64_t nv = maxvtx + 1;

	int64_t k, k1, k2, newk2;
	int64_t * vlist;


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
				const int64_t newv = ((p % 2) ? (*edges)[p / 2].second : (*edges)[p / 2].first);//get_v1_from_edge(&IJ[p / 2]) : get_v0_from_edge(&IJ[p / 2]));
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
void search_boost(
	int NBFS, int64_t nvtx_scale, int64_t* bfs_root, int64_t nedge, Edges* edges,
	int64_t maxvtx, int64_t* head, int64_t* next)
{
	//bool VERBOSE = true;
	int m, err = 0;
	int64_t k, t;
	int64_t* bfs_nedge = new int64_t[nvtx_scale];
	for (m = 0; m < NBFS; ++m) {
		int64_t *bfs_tree, max_bfsvtx;

		/* Re-allocate. Some systems may randomize the addres... */
		bfs_tree = new int64_t[nvtx_scale];
		assert(bfs_root[m] < nvtx_scale);

		if (VERBOSE) fprintf(stderr, "Running bfs %d...", m);

		err = make_bfs_tree_boost(bfs_tree, bfs_root[m], (edges), maxvtx, head, next);

		if (VERBOSE) fprintf(stderr, "done\n");

		if (err) {
			perror("make_bfs_tree failed");
			abort();
		}

		if (VERBOSE) fprintf(stderr, "Verifying bfs %d...", m);
		bfs_nedge[m] = verify_bfs_tree_boost(bfs_tree, maxvtx, bfs_root[m], (edges), nedge);
		if (VERBOSE) fprintf(stderr, "done\n");
		if (bfs_nedge[m] < 0) {
			fprintf(stderr, "bfs %d from %" PRId64 " failed verification (%" PRId64 ")\n",
				m, bfs_root[m], bfs_nedge[m]);
			abort();
		}

		free(bfs_tree);
	}
}

///////////////////////////////////////////////////////////////////////////////
int init()//hpx_main(boost::program_options::variables_map &vm)
{
    packed_edge * IJ;
    int64_t nedge = getedge();
    uint64_t scale = getscale();
    IJ = new packed_edge[nedge];

    //makeEdgeList(IJ);
    int64_t desired_nedge = nedge;
    uint64_t userseed = 22222;
    make_graph((int)scale, desired_nedge, userseed, userseed, &nedge, (packed_edge**)(&IJ));
    
	int64_t maxvtx, maxdeg;
	get_maxedge(IJ, nedge, maxvtx);
	int64_t* deg;
	int64_t* next;
	int64_t* head = new int64_t[(2 * (maxvtx + 1) + 2 * nedge)];


	deg = &head[maxvtx + 1];
	next = &deg[maxvtx + 1];

	for (int64_t k = 0; k <= maxvtx; ++k) {
		head[k] = -1;
		deg[k] = 0;
	}

	create_graph_from_edgelist(IJ, nedge, maxvtx, maxdeg, *(&head), *(&deg), *(&next));


    int64_t* bfs_root = new int64_t[getNBFS_max()];

    verify(getNBFS_max(), getnvtx_scale(), getNBFS(), nedge, IJ, bfs_root);
	search(getNBFS(), getnvtx_scale(), bfs_root, nedge, IJ, maxvtx, head, next);
	if (true)
	{
		using namespace boost;
		//graph_t g;
		Edges edgelist = Edges(nedge);
		for (int i = 0; i < nedge; ++i)
		{
			edgelist[i] = Edge(IJ[i].v0_low, IJ[i].v1_low);
			//add_edge(IJ[i].v0_low, IJ[i].v1_low, g);
		}

		int64_t* bfs_root2 = new int64_t[getNBFS_max()];
		verify_boost(getNBFS_max(), getnvtx_scale(), getNBFS(), nedge, &edgelist, bfs_root2);
		search_boost(getNBFS_max(), getnvtx_scale(), bfs_root2, nedge, &edgelist, maxvtx, head, next);
	}
    return 0;
}
int main() {
    return init();//hpx::init();
}