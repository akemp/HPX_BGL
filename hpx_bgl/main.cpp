//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "headers.hpp"

typedef std::pair < idx_t, idx_t > Edge;
typedef std::vector<Edge> Edges;

int getres(std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, std::vector<idx_t> &part,
	idx_t nparts)
{

	// Needed by parmetis
	idx_t nvtxs = xadj.size() - 1, ncon = 1;
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

typedef adjacency_list < vecS, vecS, undirectedS, vpartition > graph_t;

typedef graph_traits < graph_t >::vertices_size_type Size;

typedef
iterator_property_map < std::vector<Size>::iterator,
property_map<graph_t, vertex_index_t>::const_type >
bfstree_pm_type;

typedef hpx::lcos::local::spinlock mutex_type;
struct vgraph
{
	Size part;
	graph_t graph;
	std::vector < Size >::iterator bfstree;
	std::vector<boost::default_color_type>* colorMap;
	Size start;
	vgraph* self;
	mutex_type m;
};

struct vpartition
{
	Size part;
	vgraph* graph;
	bool visited = false;
	std::shared_ptr<vgraph> getgraph()
	{
		return std::shared_ptr<vgraph>(graph);
	}
};


template < typename TimeMap > class bfs_time_visitor :public default_bfs_visitor {
	typedef typename property_traits < TimeMap >::value_type T;
public:
	bfs_time_visitor(TimeMap tmap, T & t, Size s) :m_timemap(tmap), m_time(t), starter(s) { }

	template < typename Edge, typename Graph >
	void tree_edge(Edge e, const Graph & g)
	{
		Size src = source(e, g);
		Size targ = target(e, g);
		Size checker = get(m_timemap, targ);
		Size start = starter;
		if (checker >= num_vertices(g))
		{
			put(m_timemap, targ, src); // put source vertex in target vertex
		}
		else
		{
			runTest(checker, start, src, targ);
		}
	}
	template < typename Vertex, typename Graph >
	void finish_vertex(Vertex u, const Graph & g)
	{
		if (g[u].part != m_time)
		{
			hpx::async(&runGraph, u, g[u], starter);
		}
	}
	void runTest(Size checker, Size start, Size src, Size targ)
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
	TimeMap m_timemap;
	T & m_time;
	Size starter;
};
//runGraph(start, loc, graphs[loc], start);
int runGraph(
	Size start,
	vpartition part,
	Size begin)
{
	vgraph* g = part.graph;
	mutex_type::scoped_lock l(g->m);

	//std::atomic_store<int>(&locked, std::make_shared<int>(1));
	bfstree_pm_type parents((g->bfstree), get(vertex_index, g->graph));
	bfs_time_visitor < bfstree_pm_type >vis(parents, part.part, begin);
	typedef
		iterator_property_map < std::vector<boost::default_color_type>::iterator,
		property_map<graph_t, vertex_index_t>::const_type >
		colormap_t;
	//(*g->colorMap)[begin] = boost::default_color_type::black_color;
	colormap_t colors(g->colorMap->begin(), get(vertex_index, g->graph)); //Create a color map

	breadth_first_visit(g->graph, vertex(start, g->graph), visitor(vis).vertex_color_map(colors));
	return 0;
}


inline void runBFS(std::vector<std::vector<idx_t>> &nodes, std::vector<int>& part, int start, int partitions,
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
	std::vector<vgraph> graphs(partitions);
	vector<vector<boost::default_color_type>> map(
		partitions,
		vector<boost::default_color_type>(num_vertices(g), boost::default_color_type::white_color));
	for (int k = 0; k < partitions; ++k)
	{
		graphs[k].part = k;
		graphs[k].graph = graph_t(num_vertices(g));
		graphs[k].start = start;
		graphs[k].bfstree = (bfstree.begin());
		graphs[k].colorMap = &map[k];
		graphs[k].self = &(graphs[k]);
		vpartition part;
		part.graph = &graphs[k];
		part.part = k;
		std::fill(graphs[k].graph.m_vertices.begin(), graphs[k].graph.m_vertices.end(), part);
	}
	for (int k = 0; k < partitions; ++k)
	{
		cout << k << endl;
		for (std::size_t i = 0; i < edges.size(); ++i)
		{
			Edge e = edges[i];
			if (g[e.first].part == k && g[e.second].part == k)
			{
				add_edge(e.first, e.second, graphs[k].graph);
				//add_edge(e.second, e.first, graphs[k].graph);
			}
			else if (g[e.first].part == k)
			{
				add_edge(e.first, e.second, graphs[k].graph);
				Size part = g[e.second].part;
				graphs[k].graph[e.second].graph = graphs[part].self;
				graphs[k].graph[e.second].part = part;
			}

			else if (g[e.second].part == k)
			{
				add_edge(e.second, e.first, graphs[k].graph);
				Size part = g[e.first].part;
				graphs[k].graph[e.first].graph = graphs[part].self;
				graphs[k].graph[e.first].part = part;
			}
		}
	}

	Size loc = g[start].part;
	hpx::util::high_resolution_timer t;
	runGraph(start, graphs[loc].graph[start], start);
	hpx::finalize();
	cout << "time taken: " << t.elapsed() << "s" << " code " << endl;

	boost::random::mt19937 rng;
	boost::random::uniform_int_distribution<> randnodes(0, nodes.size() - 1);
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
int main()
{
	
	using namespace std;
	cout << "threads: " << hpx::get_os_thread_count() << endl;
	vector<idx_t> xadj;
	vector<idx_t> adjncy;

	int nnodes;
	cout << "Enter nodes: ";
	cin >> nnodes;
	int ind;
	cout << "Enter edges per node: ";
	cin >> ind;

	vector<vector<idx_t>> nodes(nnodes);
	boost::random::mt19937 rng;
	boost::random::uniform_int_distribution<> randnodes(0, nodes.size()-1);
	for (int i = 0; i < nodes.size(); ++i)
	{
		vector<idx_t> edger;
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

	if (partitions > 1)
		getres(xadj, adjncy, part, partitions);
	else
	{
		partitions = 1;
		std::cout << "Serial execution.\n";
		std::fill(part.begin(), part.end(), 0);
	}

	//vector<hpx::future<void>> futures;

	vector<Size> bfstree(nodes.size(), nodes.size());
	//for (int i = 0; i < 1; ++i)
	Size randomval = randnodes(rng);
	runBFS(nodes, part, randomval, partitions, bfstree);
	/*std::fill(part.begin(), part.end(), 0);
	std::fill(bfstree.begin(), bfstree.end(), nodes.size());
	cout << "Serial:\n";
	runBFS(nodes, part, randomval, 1, bfstree);*/
	int s;
	cin >> s;
	exit(0);
	return 0;
}