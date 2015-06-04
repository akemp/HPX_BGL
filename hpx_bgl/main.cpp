//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include "headers.hpp"

#include <queue>
#include <hpx/lcos/local/composable_guard.hpp>
#include <hpx/include/parallel_algorithm.hpp>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>

// Define the type of the graph - this specifies a bundled property for vertices

struct NodeInfo
{
	int part;
	int comp = -1;
	int dist = 999999;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, NodeInfo> Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
template < typename TimeMap > class custom_bfs_visitor : public boost::default_bfs_visitor
{
public:

	custom_bfs_visitor(TimeMap tmap, TimeMap dmap) :m_timemap(tmap), m_distmap(dmap) { }
	template < typename Edge, typename Graph >
	void tree_edge(Edge e, const Graph & g)
	{
		using namespace std;
		int target = boost::target(e, g);
		int source = boost::source(e, g);
		put(m_timemap, target, source);
		int sdist = get(m_distmap, source);
		put(m_distmap, target, sdist + 1);
	}
	TimeMap m_timemap;
	TimeMap m_distmap;
};

typedef
boost::iterator_property_map < std::vector<int>::iterator,
boost::property_map<Graph, boost::vertex_index_t>::const_type >
dtime_pm_type;
struct queueLock
{
	std::deque<int> q;
	typedef hpx::lcos::local::spinlock mutex_type;
	int front()
	{
		mutex_type::scoped_lock g(*mtx_);
		return q.front();
	}
	void push(int i)
	{
		mutex_type::scoped_lock g(*mtx_);
		q.push_back(i);
	}
	void pop()
	{
		mutex_type::scoped_lock g(*mtx_);
		q.pop_front();
	}
	bool empty()
	{
		mutex_type::scoped_lock g(*mtx_);
		return q.empty();
	}
	queueLock()
	{
		mtx_ = boost::shared_ptr<mutex_type>(new mutex_type);
	}
	boost::shared_ptr<mutex_type> mtx_;
};

void BFS_act(Graph& g, std::vector<queueLock>& q, int level, int part, bool& contains)
{
	using namespace std;
	if (q[part].empty())
	{
		return;
	}
	contains = true;
	while (!q[part].empty())
	{
		int ind = q[part].front();
		int cdist = g[ind].dist;
		if (cdist > level)
		{
			return;
		}
		q[part].pop();
		Graph::adjacency_iterator neighbourIt, neighbourEnd;
		tie(neighbourIt, neighbourEnd) = adjacent_vertices(ind, g);
		for (; neighbourIt != neighbourEnd; ++neighbourIt)
		{
			int test = *neighbourIt;
			if (g[test].dist > cdist + 1)
			{
				g[test].dist = cdist + 1;
				g[test].comp = ind;
				int p = g[test].part;
				q[p].push(test);
			}
		}
	}
	return;
}

void BFS_Search(Graph& g, int start, int parts)
{
	using namespace boost;
	using namespace hpx;
	g[start].dist = 0;
	std::vector<queueLock> q(parts);
	int stpart = g[start].part;
	q[stpart].push(start);
	int level = 0;
	bool cont = true;
	
	while (cont)
	{
		using namespace std;
		cont = false;
		{
			typedef boost::counting_iterator<std::size_t> iterator;
			hpx::parallel::for_each(hpx::parallel::par, iterator(0), iterator(q.size()),
				[&g,&q,level,&cont](std::size_t i){
				BFS_act(g, q, level, i, cont); 
			});
		}
		++level;
	}
	
}

int main()
{

	using namespace std;
	cout << "Running on " << hpx::get_os_thread_count() << " threads.\n";
    uint64_t nnodes;
    uint64_t scale;
	cout << "Enter scale (nodes=2^input): ";
	cin >> scale;
    uint64_t edgefactor;
	cout << "Enter edgefactor: ";
    cin >> edgefactor;

	int searches = 1;
    nnodes = (uint64_t)(1) << scale;
	//


	//Edges edges;
	vector<vector<int>> nodes(
		nnodes
		);
    uint64_t ind = (uint64_t)(edgefactor) << scale;
	vector<int> pedges(ind);
	boost::random::mt19937 rng;
	boost::random::uniform_int_distribution<> randnodes(0, nnodes-1);
	for (int i = 0; i < pedges.size(); ++i)
		pedges[i] = randnodes(rng);
	cout << "Random range generated. Making edgelist.\n";
	{
		parallel_edge_gen(pedges.begin(),nodes, ind);
	}

	cout << "Edgelist generated. Creating partitions.\n";
	vector<int> xadj, adjncy;
	toCSR(nodes, xadj, adjncy);
	vector<int> parts(nodes.size());
	int nparts = 8;
	get_parts(xadj, adjncy, parts, nparts);
	cout << "Running tests.\n";

	vector<int> starts(searches);

	cout << "Setting up serial values for testing.\n";
	cout << "Data structures set. Running serial search.\n";

	//if (false)

	std::vector < int > serialcomp(nnodes, -1);
	std::vector < int > sdists(nnodes, 99999);
	{
		using namespace boost;
		Graph G(nnodes);
		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int j = 0; j < nodes[i].size(); ++j)
				boost::add_edge(i, nodes[i][j], G);
		}


		dtime_pm_type dtime_pm(serialcomp.begin(), get(vertex_index, G));
		typedef
			iterator_property_map < std::vector<int>::iterator,
			property_map<Graph, vertex_index_t>::const_type >
			dtime_pm_type;
		sdists[0] = 0;

		dtime_pm_type dists_pm(sdists.begin(), get(vertex_index, G));

		custom_bfs_visitor < dtime_pm_type >vis(dtime_pm, dists_pm);

		hpx::util::high_resolution_timer t;

		breadth_first_search(G, boost::vertex(0, G), visitor(vis));

		cout << t.elapsed() << "s\n" << endl;
	}

	cout << "Serial search completed. Running parallel search.\n";
	Graph G(nnodes);
	{
		using namespace boost;
		using namespace std;

		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int j = 0; j < nodes[i].size(); ++j)
				boost::add_edge(i, nodes[i][j], G);
		}
		for (int i = 0; i < nnodes; ++i)
		{
			G[i].part = parts[i];
		}
		hpx::util::high_resolution_timer t;

		BFS_Search(G, 0, nparts);

		cout << t.elapsed() << "s\n" << endl;
		
	}
	hpx::wait_all();
	cout << "Parallel search completed. Checking accuracy.\n";
	for (int i = 0; i < serialcomp.size(); ++i)
	{
		if (sdists[i] != G[i].dist)
		{
			std::cout << "ERROR! Serial distance of ";
			std::cout << sdists[i] << " != parallel distance of " << G[i].dist << std::endl;
			break;
		}
	}
	system("Pause");
	return 0;
}