/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <fstream>

#include <assert.h>

//#include <alloca.h> /* Portable enough... */
#include <fcntl.h>
//#include <unistd.h>

#if !defined(__MTA__)
//#include <getopt.h>
#endif


#include "graph500.h"
#include "rmat.h"
#include "kronecker.h"
#include "verify.h"
#include "prng.h"
#include "xalloc.h"
#include "options.h"
#include "generator/splittable_mrg.h"
#include "generator/make_graph.h"

static int64_t bfs_root[NBFS_max];

int getedge()
{
    return nvtx_scale * edgefactor;
}

int
makeEdgeList(packed_edge * IJ)
{
    int64_t desired_nedge = nvtx_scale * edgefactor;
  int * restrict has_adj;
  int fd;
  if (sizeof (int64_t) < 8) {
    fprintf (stderr, "No 64-bit support.\n");
    return 1;
  }


  init_random ();

  /* Catch a few possible overflows. */
  assert (desired_nedge >= nvtx_scale);
  assert (desired_nedge >= edgefactor);


  //if (1) fprintf (stderr, "Generating edge list...");
  int64_t nedge = desired_nedge;
  rmat_edgelist(IJ, nedge, SCALE, A, B, C);  
     {
      has_adj = new int[nvtx_scale];
      OMP("omp parallel") {
          OMP("omp for")
              for (int64_t k = 0; k < nvtx_scale; ++k)
                  has_adj[k] = 0;
          MTA("mta assert nodep") OMP("omp for")
              for (int64_t k = 0; k < nedge; ++k) {
              const int64_t i = get_v0_from_edge(&IJ[k]);
              const int64_t j = get_v1_from_edge(&IJ[k]);
              if (i != j)
                  has_adj[i] = has_adj[j] = 1;
              }
      }

      /* Sample from {0, ..., nvtx_scale-1} without replacement. */
      int m = 0;
      int64_t t = 0;
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
  }

    free(has_adj);
  return nedge;
}