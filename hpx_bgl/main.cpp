//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <hpx/hpx_init.hpp>
#include <hpx/include/runtime.hpp>
#include <hpx/include/thread_executors.hpp>
#include <hpx/util/high_resolution_timer.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include "../metis/include/metis.h"

#include <boost/random.hpp>
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

typedef struct vpartition;

typedef adjacency_list < vecS, vecS, directedS, vpartition > graph_t;

typedef graph_traits < graph_t >::vertices_size_type Size;

typedef
iterator_property_map < std::vector<Size>::iterator,
property_map<graph_t, vertex_index_t>::const_type >
bfstree_pm_type;

void runGraph(
	Size child,
	Size part,
	graph_t* graph,
	std::vector < Size >* bfstree,
	std::vector<boost::default_color_type>* colorMap,
	int* locked,
	Size start);
template < typename TimeMap > class bfs_time_visitor :public default_bfs_visitor {
	typedef typename property_traits < TimeMap >::value_type T;
public:
	bfs_time_visitor(TimeMap tmap, T & t, Size s) :m_timemap(tmap), m_time(t), starter(s) { }

	template < typename Edge, typename Graph >
	void examine_edge(Edge e, const Graph & g) 
	{
		Size src = source(e, g);
		Size targ = target(e, g);
		Size checker = get(m_timemap, targ);
		Size start = g[src].start;
		if (checker >= num_vertices(g))
		{
			put(m_timemap, targ, src); // put source vertex in target vertex
		}
		else
		{
			int count1 = 0;
			int count2 = 0;
			Size tester = checker;
			while (tester != start)
			{
				tester = get(m_timemap, tester);
				++count1;
			}
			tester = src;
			while (tester != start)
			{
				tester = get(m_timemap, tester);
				++count2;
			}
			if (count1 > count2)
				put(m_timemap, targ, src); // put source vertex in target vertex
		}
	}
	template < typename Vertex, typename Graph >
	void finish_vertex(Vertex u, const Graph & g)
	{
		Size targ = u;
		if (g[u].part != m_time)
		{
			//hpx::async(&runGraph,targ, g[targ].part, g[targ].graph, g[targ].bfstree, g[targ].colorMap, g[targ].locked, g[targ].start);
			runGraph( targ, g[targ].part, g[targ].graph, g[targ].bfstree, g[targ].colorMap, g[targ].locked, g[targ].start);
		}
	}
	TimeMap m_timemap;
	T & m_time;
	Size starter;
	std::vector<int> touse;
	bool used = true;
};

void runGraph(
	Size child,
	Size part,
	graph_t* graph,
	std::vector < Size >* bfstree,
	std::vector<boost::default_color_type>* colorMap,
	int* locked,
	Size start)
{
	*locked = 1;
	bfstree_pm_type parents(bfstree->begin(), get(vertex_index, *graph));
	bfs_time_visitor < bfstree_pm_type >vis(parents, part, start);
	typedef
		iterator_property_map < std::vector<boost::default_color_type>::iterator,
		property_map<graph_t, vertex_index_t>::const_type >
		colormap_t;

	colormap_t colors(colorMap->begin(), get(vertex_index, *graph)); //Create a color map

	breadth_first_visit(*graph, vertex(child, *graph), visitor(vis).vertex_color_map(colors));
	*locked = 0;
	return;
}

struct vpartition
{
	Size part;
	bool present;
	graph_t* graph;
	std::vector < Size >* bfstree;
	std::vector<boost::default_color_type>* colorMap;
	Size start;
	int* locked;
};
void runPartition(graph_t g, std::vector < Size >& bfstree, std::vector<Edge>& edges,int start, int partitions)
{
	using namespace std;
	std::vector<boost::default_color_type> colorMap(bfstree.size());
	vector<graph_t> graphs(partitions);
	vector<int> locks(partitions, 0);
	for (int k = 0; k < partitions; ++k)
	{
		graphs[k] = g;
		locks[k] = false;
		vpartition part;
		part.part = k;
		part.start = start;
		part.present = 0;
		part.graph = &graphs[k];
		part.bfstree = &bfstree;
		part.colorMap = &colorMap;
		part.locked =  &(locks[k]);
		std::fill(graphs[k].m_vertices.begin(), graphs[k].m_vertices.end(), part);
	}
	for (int k = 0; k < partitions; ++k)
	{
		for (std::size_t i = 0; i < edges.size(); ++i)
		{
			Edge e = edges[i];
			if (g[e.first].part == k && g[e.second].part == k)
			{
				add_edge(e.first, e.second, graphs[k]);
				add_edge(e.second, e.first, graphs[k]);
			}
			else if (g[e.first].part == k)
			{
				add_edge(e.first, e.second, graphs[k]);
				Size part = g[e.second].part;
				graphs[k][e.second].graph = &(graphs[part]);
				graphs[k][e.second].part = part;
			}

			else if (g[e.second].part == k)
			{
				add_edge(e.second, e.first, graphs[k]);
				Size part = g[e.first].part;
				graphs[k][e.first].graph = &(graphs[part]);
				graphs[k][e.first].part = part;
			}
		}
	}
	
	Size loc = g[start].part;
	runGraph(start, loc, graphs[loc][start].graph, graphs[loc][start].bfstree, graphs[loc][start].colorMap, graphs[loc][start].locked, start);
	hpx::wait_all();
	return;
}

void runBFS(std::vector<std::vector<idx_t>> &nodes, std::vector<int>& part, int start, int partitions,
	std::vector < Size >& bfstree)
{
	using namespace std;

	vector<Edge> edges;
	createEdges(nodes, edges);

	// a vector to hold the parent property for each vertex


	int size = nodes.size();
	graph_t g(size);
	for (int i = 0; i < nodes.size(); ++i)
	{
		g[i].part = part[i];
	}

	bfstree[start] = start;
	runPartition( g, bfstree, edges, start, partitions);

	boost::random::mt19937 rng;
	boost::random::uniform_int_distribution<> randnodes(0, nodes.size() - 1);
	/**/
	for (int i = 0; i < 32; ++i)
	{
		int sample = randnodes(rng);
		cout << sample << " ";
		while (sample != start)
		{
			if (sample >= nodes.size())
			{
				cout << "\nERROR!\n";
				break;
			}
			sample = bfstree[sample];
			cout << sample << " ";
		}
		cout << endl;
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////
int hpx_main()
{
	using namespace std;
	int result;
	cout << "threads: " << hpx::get_os_thread_count() << endl;
	vector<idx_t> xadj;
	vector<idx_t> adjncy;

	vector<vector<idx_t>> nodes(8192);
	boost::random::mt19937 rng;
	boost::random::uniform_int_distribution<> randnodes(0, nodes.size()-1);
	for (int i = 0; i < nodes.size(); ++i)
	{
		vector<idx_t> edger;
		int ind = 200;
		for (int j = 0; j < ind; ++j)
		{
			int spot = randnodes(rng);
			if (spot != i && std::find(nodes[spot].begin(), nodes[spot].end(), i) == nodes[spot].end())
				edger.push_back(spot);
		}
		nodes[i] = edger;
	}
	toCSR(nodes, xadj, adjncy);


	vector<idx_t> part(xadj.size()-1);
	
	int partitions = 30;
	cout << endl << "Enter partitions: ";
	cin >> partitions;

	int res = getres(xadj, adjncy, part, partitions);

	hpx::util::high_resolution_timer t;
	vector<hpx::future<void>> futures;

	vector<Size> bfstree(nodes.size(), nodes.size());
	//for (int i = 0; i < 1; ++i)
	Size randomval = randnodes(rng);
	runBFS(nodes, part, randomval, partitions, bfstree);
	cout << "time taken: " << t.elapsed() << "s" << endl;
	std::fill(part.begin(), part.end(), 0);
	std::fill(bfstree.begin(), bfstree.end(), nodes.size());
	cout << "Serial:\n";
	runBFS(nodes, part, randomval, partitions, bfstree);
	return hpx::finalize();
}

int main(int argc, char* argv[]) 
{
	return hpx::init(argc, argv);
}