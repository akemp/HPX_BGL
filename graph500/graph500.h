/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#if !defined(GRAPH500_HEADER_)
#define GRAPH500_HEADER_

#define NAME "Graph500 sequential list"
#define VERSION 0

#include "generator/graph_generator.h"
#include "verify.h"

/** Pass the edge list to an external graph creation routine. */
int create_graph_from_edgelist(struct packed_edge *IJ, int64_t nedge, int64_t &maxvtx, int64_t& maxdeg, int64_t* head);

/** Create the BFS tree from a given source vertex. */
int make_bfs_tree (int64_t *bfs_tree_out, int64_t *max_vtx_out,
	int64_t srcvtx, packed_edge * IJ, int64_t maxvtx, int64_t* head);

void get_maxedge(struct packed_edge *IJ_in, int64_t nedge, int64_t &maxvtx, int64_t& maxdeg);
/** Clean up. */
void destroy_graph (void);

void
run_bfs(packed_edge * IJ, int NBFS_max, int64_t nvtx_scale, int NBFS, int64_t nedge);

#endif /* GRAPH500_HEADER_ */
