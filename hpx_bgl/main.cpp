//  Copyright (c) 2011 Matthew Anderson
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <hpx/hpx_init.hpp>
#include <hpx/include/runtime.hpp>
#include <hpx/include/thread_executors.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include "../metis/include/metis.h"

#include <boost/shared_array.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>

//int64_t *deg, *next;
typedef std::pair < idx_t, idx_t > Edge;
typedef std::vector<Edge> Edges;

int getres(std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, std::vector<idx_t> &part,
	idx_t nparts)
{

	// Needed by parmetis
	idx_t nvtxs = xadj.size()-1, ncon = 1;
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

void createEdges(const std::vector<std::vector<idx_t>>& nodes, std::vector<Edge>& edges)
{

	for (int i = 0; i < nodes.size(); ++i)
	{
		for (int j = 0; j < nodes[i].size(); ++j)
		{
			edges.push_back(Edge(i, nodes[i][j]));
		}
	}
	return;
}

using namespace boost;

struct partition
{
	int part;
};

typedef adjacency_list < vecS, vecS, undirectedS, partition > graph_t;


template < typename TimeMap > class bfs_time_visitor :public default_bfs_visitor {
	typedef typename property_traits < TimeMap >::value_type T;
public:
	bfs_time_visitor(TimeMap tmap, T & t) :m_timemap(tmap), m_time(t) { }
	/*template < typename Vertex, typename Graph >
	void discover_vertex(Vertex u, const Graph & g) const
	{
		int val = get(m_timemap, u);
		if (val != m_time)
		{
			std::cout << "Changing indices to " << val << " from " << m_time << std::endl;
		}
	}*/
	template < typename Edge, typename Graph >
	void gray_target(Edge e, const Graph & g) const
	{
		/*int src = get(m_timemap, source(e, g));
		int targ = get(m_timemap, target(e, g));
		if (src != targ)
		{

		}*/
		put(m_timemap, source(e, g), target(e, g));
	}
	TimeMap m_timemap;
	T & m_time;
};

void runBFS(std::vector<std::vector<idx_t>> &nodes, std::vector<Edge>& edges, int start)
{

	{
		graph_t g(nodes.size());
		using namespace boost;
		using namespace std;
		// Select the graph type we wish to use
		// Create the graph object
		const int n_edges = edges.size();//sizeof(edge_array) / sizeof(E);
		for (std::size_t j = 0; j < n_edges; ++j)
			add_edge(edges[j].first, edges[j].second, g);

		// Typedefs
		typedef graph_traits < graph_t >::vertices_size_type Size;

		// a vector to hold the discover time property for each vertex
		std::vector < Size > bfstree(num_vertices(g));
		//for (int i = 0; i < bfstree.size(); ++i)
		//	bfstree[i] = part[i];

		typedef
			iterator_property_map < std::vector<Size>::iterator,
			property_map<graph_t, vertex_index_t>::const_type >
			bfstree_pm_type;
		bfstree_pm_type parents(bfstree.begin(), get(vertex_index, g));

/*		for (int i = 0; i < num_vertices(g); ++i)
		{
			g[i].part = part[i];
		}
*/
		Size time = 0;
		bfs_time_visitor < bfstree_pm_type >vis(parents, time);
		breadth_first_search(g, vertex(start, g), visitor(vis));
		//for (int i = 0; i < bfstree.size(); ++i)
		//	cout << bfstree[i] << " ";

	}
}

///////////////////////////////////////////////////////////////////////////////
int hpx_main()
{
	using namespace std;
	int result;
	cout << hpx::get_os_thread_count() << endl;
	//from librelocus

	vector<idx_t> xadj;
	vector<idx_t> adjncy;

	vector<vector<idx_t>> nodes(4096);
	for (int i = 0; i < nodes.size(); ++i)
	{
		vector<idx_t> edger;
		int ind = rand() % 24;
		for (int j = 0; j < ind; ++j)
		{
			edger.push_back(rand() % (nodes.size()-1));
		}
		nodes[i] = edger;
	}
	toCSR(nodes, xadj, adjncy);


	vector<idx_t> part(xadj.size());
	int res = getres(xadj, adjncy, part, 16);



	vector<Edge> edges;
	createEdges(nodes, edges);

	for (int i = 0; i < 4; ++i)
		hpx::async(hpx::util::bind(&runBFS, nodes, edges, rand()%nodes.size()));
	hpx::wait_all();
	return hpx::finalize();
}
int main(int argc, char* argv[]) {
	return hpx::init(argc, argv);
}