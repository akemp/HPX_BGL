/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include <assert.h>

//#include <alloca.h> /* Portable enough... */
#include <fcntl.h>
/* getopt should be in unistd.h */
//#include <unistd.h>

#if !defined(__MTA__)
//#include <getopt.h>
#endif

#include "graph500.h"
#include "rmat.h"
#include "kronecker.h"
#include "verify.h"
#include "prng.h"
//#include "timer.h"
#include "xalloc.h"
//#include "options.h"
#include "generator/splittable_mrg.h"
#include "generator/graph_generator.h"
#include "generator/make_graph.h"

//static const struct packed_edge * restrict IJ;



/*
void
run_bfs(packed_edge * IJ, int NBFS_max, int64_t nvtx_scale, int NBFS, int64_t nedge)
{

    int64_t* bfs_root = new int64_t[NBFS_max];

    static double generation_time;
    static double construction_time;
    double* bfs_time = new double[NBFS_max];
    int64_t* bfs_nedge = new int64_t[NBFS_max];

    int * has_adj;
    int m, err;
    int64_t k, t;
    
    if (VERBOSE) fprintf (stderr, "Creating graph...");
        err = create_graph_from_edgelist (IJ, nedge);
    if (VERBOSE) fprintf (stderr, "done.\n");
    if (err) {
        fprintf (stderr, "Failure creating graph.\n");
        exit (EXIT_FAILURE);
    }

    if (1) { //if (!rootname) {
        has_adj = new int[nvtx_scale];
        OMP("omp parallel") {
            OMP("omp for")
                for (k = 0; k < nvtx_scale; ++k)
                    has_adj[k] = 0;
            MTA("mta assert nodep") OMP("omp for")
                for (k = 0; k < nedge; ++k) {
                const int64_t i = get_v0_from_edge(&IJ[k]);
                const int64_t j = get_v1_from_edge(&IJ[k]);
                if (i != j)
                    has_adj[i] = has_adj[j] = 1;
                }
        }

        m = 0;
        t = 0;
        while (m < NBFS && t < nvtx_scale) {
            double R = mrg_get_double_orig (prng_state);
            if (!has_adj[t] || (nvtx_scale - t)*R > NBFS - m) ++t;
            else bfs_root[m++] = t++;
        }
        if (t >= nvtx_scale && m < NBFS) {
            if (m > 0) {
                fprintf (stderr, "Cannot find %d sample roots of non-self degree > 0, using %d.\n",
                    NBFS, m);
                NBFS = m;
            } else {
                fprintf (stderr, "Cannot find any sample roots of non-self degree > 0.\n");
                exit (EXIT_FAILURE);
            }
        }

        free (has_adj);
    }

    for (m = 0; m < NBFS; ++m) {
        int64_t *bfs_tree, max_bfsvtx;

        bfs_tree = new int64_t[nvtx_scale];
        assert (bfs_root[m] < nvtx_scale);

        if (VERBOSE) fprintf (stderr, "Running bfs %d...", m);
            err = make_bfs_tree(bfs_tree, &max_bfsvtx, bfs_root[m], &(*IJ));//TIME(bfs_time[m], err = make_bfs_tree (bfs_tree, &max_bfsvtx, bfs_root[m]));
        if (VERBOSE) fprintf (stderr, "done\n");

        if (err) {
            perror ("make_bfs_tree failed");
            abort ();
        }

        if (VERBOSE) fprintf (stderr, "Verifying bfs %d...", m);
        bfs_nedge[m] = verify_bfs_tree (bfs_tree, max_bfsvtx, bfs_root[m], IJ, nedge);
        if (VERBOSE) fprintf (stderr, "done\n");
        if (bfs_nedge[m] < 0) {
            fprintf (stderr, "bfs %d from %" PRId64 " failed verification (%" PRId64 ")\n",
                m, bfs_root[m], bfs_nedge[m]);
            abort ();
        }

        free (bfs_tree);
    }

    destroy_graph ();
}
*/