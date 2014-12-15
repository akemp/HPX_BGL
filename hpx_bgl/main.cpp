//  Copyright (c) 2011 Matthew Anderson
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include "../metis/include/metis.h"
//#include <hpx/hpx.hpp>
//#include <hpx/hpx_init.hpp>

#include <hpx/util/high_resolution_timer.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/range/irange.hpp>


typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS > graph_t;

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
template < typename TimeMap > class bfs_time_visitor :public default_bfs_visitor {
	typedef typename property_traits < TimeMap >::value_type T;
public:
	bfs_time_visitor(TimeMap tmap, T & t) :m_timemap(tmap), m_time(t) { }
	template < typename Vertex, typename Graph >
	void discover_vertex(Vertex u, const Graph & g) const
	{
		put(m_timemap, u, m_time++);
	}
	TimeMap m_timemap;
	T & m_time;
};

///////////////////////////////////////////////////////////////////////////////
int init()//hpx_main(boost::program_options::variables_map &vm)
{
	using namespace std;
	int result;
	//from librelocus

	vector<idx_t> xadj;
	vector<idx_t> adjncy;

	vector<vector<idx_t>> nodes(1024);
	for (int i = 0; i < nodes.size(); ++i)
	{
		vector<idx_t> edger;
		int ind = rand() % 12;
		for (int j = 0; j < ind; ++j)
		{
			edger.push_back(rand() % (nodes.size()-1));
		}
		nodes[i] = edger;
	}
	toCSR(nodes, xadj, adjncy);

	vector<idx_t> part(xadj.size());
	int res = getres(xadj, adjncy, part, 16);
	/*cout << res << endl;
	cout << "npart is " << '\n';
	for (int i = 0; i<part.size(); i++){
		cout << part[i] << " ";
	}
	cout << '\n';*/
	vector<Edge> edges;
	createEdges(nodes, edges);
	graph_t g(nodes.size());


	cout << "Metis returned successfully." << '\n';

	{
		using namespace boost;
		// Select the graph type we wish to use
		typedef adjacency_list < vecS, vecS, undirectedS > graph_t;
		// Set up the vertex IDs and names
		const char *name = "rstuvwxy";
		// Specify the edges in the graph
		/*typedef std::pair < int, int >E;
		E edge_array[] = { E(r, s), E(r, v), E(s, w), E(w, r), E(w, t),
			E(w, x), E(x, t), E(t, u), E(x, y), E(u, y)
		};*/
		// Create the graph object
		const int n_edges = edges.size();//sizeof(edge_array) / sizeof(E);
		// VC++ has trouble with the edge iterator constructor
		graph_t g(nodes.size());
		for (std::size_t j = 0; j < n_edges; ++j)
			add_edge(edges[j].first, edges[j].second, g);

		// Typedefs
		typedef graph_traits < graph_t >::vertices_size_type Size;

		// a vector to hold the discover time property for each vertex
		std::vector < Size > dtime(num_vertices(g));
		typedef
			iterator_property_map < std::vector<Size>::iterator,
			property_map<graph_t, vertex_index_t>::const_type >
			dtime_pm_type;
		dtime_pm_type dtime_pm(dtime.begin(), get(vertex_index, g));

		Size time = 0;
		bfs_time_visitor < dtime_pm_type >vis(dtime_pm, time);
		breadth_first_search(g, vertex(3, g), visitor(vis));

		// Use std::sort to order the vertices by their discover time
		std::vector<graph_traits<graph_t>::vertices_size_type > discover_order(nodes.size());
		integer_range < int >range(0, nodes.size());
		std::copy(range.begin(), range.end(), discover_order.begin());
		std::sort(discover_order.begin(), discover_order.end(),
			indirect_cmp < dtime_pm_type, std::less < Size > >(dtime_pm));

		std::cout << "order of discovery: ";
		for (int i = 0; i < nodes.size(); ++i)
			std::cout << discover_order[i] << " ";
		std::cout << std::endl;
	}

	return 0;
}
int main() {
    return init();//hpx::init();
}