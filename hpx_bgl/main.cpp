//  Copyright (c) 2011 Matthew Anderson
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
//#include <hpx/hpx.hpp>
// #include <hpx/hpx_init.hpp>

#include <hpx/util/high_resolution_timer.hpp>

#include "../graph500/make-edgelist.h"
#include "../graph500/graph500.h"
#include "../graph500/generator/make_graph.h"

#include "../graph500/prng.h"
#include "../graph500/generator/splittable_mrg.h"

void verify(int NBFS_max, int64_t nvtx_scale, int NBFS, int64_t nedge,
    packed_edge* IJ, int64_t* bfs_root)
{

    int * has_adj;
    int m, err;
    int64_t k, t;

    has_adj = new int[nvtx_scale];
    {
            for (k = 0; k < nvtx_scale; ++k)
                has_adj[k] = 0;
            for (k = 0; k < nedge; ++k) {
                const int64_t i = get_v0_from_edge(&IJ[k]);
                const int64_t j = get_v1_from_edge(&IJ[k]);
                if (i != j)
                    has_adj[i] = has_adj[j] = 1;
            }
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

void search(int NBFS, int64_t nvtx_scale, int64_t* bfs_root, int64_t nedge, packed_edge* IJ)
{
    bool VERBOSE = true;
    int m, err;
    int64_t k, t;
    int64_t* bfs_nedge = new int64_t[nvtx_scale];
    for (m = 0; m < NBFS; ++m) {
        int64_t *bfs_tree, max_bfsvtx;

        /* Re-allocate. Some systems may randomize the addres... */
        bfs_tree = new int64_t[nvtx_scale];
        assert(bfs_root[m] < nvtx_scale);

        if (VERBOSE) fprintf(stderr, "Running bfs %d...", m);
        
        err = make_bfs_tree(bfs_tree, &max_bfsvtx, bfs_root[m]);//TIME(bfs_time[m], err = make_bfs_tree (bfs_tree, &max_bfsvtx, bfs_root[m]));
        
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

///////////////////////////////////////////////////////////////////////////////
int init()//hpx_main(boost::program_options::variables_map &vm)
{
    packed_edge * IJ;
    int64_t nedge = getedge();
    uint64_t scale = getscale();
    IJ = new packed_edge[nedge];

    //makeEdgeList(IJ);
    int64_t desired_nedge = nedge;
    uint64_t userseed = 1000;
    make_graph((int)scale, desired_nedge, userseed, userseed, &nedge, (packed_edge**)(&IJ));
    
    create_graph_from_edgelist(IJ, nedge);

    int64_t* bfs_root = new int64_t[getNBFS_max()];

    verify(getNBFS_max(), getnvtx_scale(), getNBFS(), nedge, IJ, bfs_root);
    search(getNBFS(), getnvtx_scale(), bfs_root, nedge, IJ);
    /*std::ofstream fout("out.txt");
    for (int i = 0; i < nedge; ++i)
    {
        fout << IJ[i].v0_low << " " << IJ[i].v1_low << std::endl;;
    }
    fout.close();
    */
    //run_bfs(IJ, getNBFS_max(), getnvtx_scale(), getNBFS(), nedge);
    //hpx::finalize();
    return 0;
}
int main() {
    return init();//hpx::init();
}