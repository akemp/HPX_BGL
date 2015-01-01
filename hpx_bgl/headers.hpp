#include <hpx/hpx_main.hpp>
#include <hpx/include/runtime.hpp>
#include <hpx/include/thread_executors.hpp>
#include <hpx/util/high_resolution_timer.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/iostreams.hpp>


#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <mutex>
#include <thread>         // std::thread
#include "../metis/include/metis.h"

#include <boost/random.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>


typedef std::pair < uint32_t, uint32_t > Edge;
typedef std::vector<Edge> Edges;
using namespace boost;
using namespace std;
