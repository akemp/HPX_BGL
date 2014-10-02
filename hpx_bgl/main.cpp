//  Copyright (c) 2011 Matthew Anderson
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>

#include <hpx/util/high_resolution_timer.hpp>


#include "../graph500/make-edgelist.h"
///////////////////////////////////////////////////////////////////////////////
int hpx_main(boost::program_options::variables_map &vm)
{
    packed_edge * IJ;
    makeEdgeList(IJ);
    std::ofstream fout("out.txt");
    for (int i = 0; i < nedge; ++i)
    {
        fout << IJ[i].v0_low << " " << IJ[i].v1_low << std::endl;;
    }
    fout.close();
    hpx::finalize();
    return 0;
}
int main() {
    return hpx::init();
}