#ifndef HEADERS_BGL
#define HEADERS_BGL

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
#include "../metis/include/metis.h"

#include <boost/random.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include "high_resolution_timer.hpp"


typedef std::pair < int, int > Edge;
typedef std::vector<Edge> Edges;

struct multi_name_t {
	typedef boost::vertex_property_tag kind;
};

typedef boost::property<multi_name_t, std::vector<int> > MultiColor; //parent, color, partition, distance
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	MultiColor> Graph;

int get_parts(std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, std::vector<idx_t> &part,
	idx_t nparts);
void toCSR(const std::vector<std::vector<idx_t>>& nodes, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy);

using namespace boost;
using namespace std;
struct GraphComponent
{
	static void parreset(Graph* g, int start, int size, int toggled)
	{
		property_map<Graph, multi_name_t>::type
			name = get(multi_name_t(), *g);
		for (int i = start; i < start + size; ++i)
		{
			name[i] = vector<int>(toggled, -1);
		}
	}

	void reset(int toggled)
	{
		int adder = grainsize;

		vector<hpx::thread> futs;
		Graph* ptr = &g;

		for (int i = 0; i < num_vertices(g); i += adder)
		{
			if (i + adder < num_vertices(g))
			{
				futs.push_back(hpx::thread(&parreset, ptr, i, adder, toggled));
			}
			else
			{
				futs.push_back(hpx::thread(&parreset, ptr, i, num_vertices(g) - i, toggled));
			}
		}
		for (int i = 0; i < futs.size(); ++i)
			futs[i].join();

	}
	void set(vector<vector<int>> nodes, int size, int edge, int starts)
	{
		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int j = 0; j < nodes[i].size(); ++j)
				add_edge(i, nodes[i][j], g);
		}
		grainsize = size;
		edgefact = edge;

		name = get(multi_name_t(), g);
		multireset(starts);
	}
	int getval(int i, int index)
	{
		return name[i][index];
	}
	void multireset(int starts)
	{
		for (int i = 0; i < num_vertices(g); ++i)
		{
			name[i] = vector<int>(starts, -1);
		}
	}

	void setmulti(vector<vector<int>> nodes, int size, int edge, int starts)
	{
		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int j = 0; j < nodes[i].size(); ++j)
				add_edge(i, nodes[i][j], g);
		}
		grainsize = size;
		edgefact = edge;

		name = get(multi_name_t(), g);
		multireset(starts);
	}

	void multival(vector<int> starts, bool sequential)
	{
		if (false)
		{
			vector<hpx::thread> futures;
			int i = 0;
			int adder = 1;
			for (vector<int>::iterator it = starts.begin(); it < starts.end(); it += adder)
			{
				int last = i;
				i += adder;
				if (i < starts.size())
				{
					futures.push_back(hpx::thread(&runMp, it, last, adder, this));
				}
				else
				{
					futures.push_back(hpx::thread(&runMp, it, last, starts.size() - last, this));
					break;
				}
			}
			for (int i = 0; i < futures.size(); ++i)
				futures[i].join();
		}
		else
		{
			for (int i = 0; i < starts.size(); ++i)
			{
				mpbfs(starts[i], i);
			}
		}
	}
	static void runMp(vector<int>::iterator loc, int index, int size, GraphComponent* gc)
	{
		for (int i = 0; i < size; ++i)
		{
			gc->mpbfs(*loc, index + i);
			++loc;
		}
	}
	static vector<int> process_layor_multi(int loc, vector<int> in_bag, Graph* g)
	{
		property_map < Graph, vertex_index_t >::type
			index_map = get(vertex_index, *g);
		property_map<Graph, multi_name_t>::type
			name = get(multi_name_t(), *g);
		vector<int> out_bag;
		int count = 0;
		for (int i = 0; i < in_bag.size(); ++i)
		{
			int val = in_bag[i];
			graph_traits < Graph >::adjacency_iterator ai, a_end;

			for (boost::tie(ai, a_end) = adjacent_vertices(val, *g); ai != a_end; ++ai)
			{
				int ind = get(index_map, *ai);
				if (name[ind][loc] >= 0)
					continue;
				name[ind][loc] = val;
				out_bag.push_back(ind);
			}
		}
		return out_bag;
	}
	void mpbfs(int index, int loc)
	{
		name[index][loc] = index;
		vector<int> v;
		int dist = 0;
		v.push_back(index);
		Graph* ptr = &g;
		while (!v.empty())
		{
			vector<hpx::future<vector<int>>> futures;
			futures.reserve(v.size() / grainsize + 1);
			{
				int i = 0;
				for (vector<int>::iterator it = v.begin(); it < v.end(); it += grainsize)
				{
					int last = i;
					i += grainsize;
					if (i < v.size())
						futures.push_back(hpx::async(std::bind(&process_layor_multi, loc, vector<int>(it, it + grainsize), ptr)));
					else
					{
						futures.push_back(hpx::async(std::bind(&process_layor_multi, loc, vector<int>(it, it + (v.size() - last)), ptr)));
						break;
					}
				}
			}
			vector<int> children;
			for (int i = 0; i < futures.size(); ++i)
			{
				vector<int> future = futures[i].get();
				children.insert(children.end(), future.begin(), future.end());
			}
			v = children;
		}
	}
	int getmultival(int index, int i)
	{
		return name[index][i];
	}

	void bfs_search_act(vector<int> starts)
	{

		for (int i = 0; i < starts.size(); ++i)
		{
			bfs_search(starts[i], i);
		}
	}

	void bfs_search(int index, int loc)
	{
		property_map < Graph, vertex_index_t >::type
			index_map = get(vertex_index, g);
		//pennants = std::vector <int>(num_vertices(g), -1);
		name[index][loc] = index;
		vector<int> q;
		q.reserve(num_vertices(g));
		q.push_back(index);
		int dist = 0;
		int spot = 0;
		while (spot < q.size())
		{
			int ind = q[spot];
			++spot;
			int parent = ind;
			graph_traits < Graph >::adjacency_iterator ai, a_end;

			for (boost::tie(ai, a_end) = adjacent_vertices(ind, g); ai != a_end; ++ai)
			{
				int ind = get(index_map, *ai);
				if (name[ind][loc] < 0)
				{
					name[ind][loc] = parent;
					q.push_back(ind);
				}
			}

		}
	};

	int getnum()
	{
		return num_vertices(g);
	}

	property_map<Graph, multi_name_t>::type
		name;
	Graph g;
	int grainsize = 9999;
	int edgefact = 16;
};

struct graph_manager
{
	GraphComponent comp;

	void pbfs_search(vector<int >index)
	{
		comp.multival(index, true);
	}
	void set(vector<vector<int>>& edges, int grainsize, int edgefact, int starts)
	{
		comp.set(edges, grainsize, edgefact, starts);
	}
	int getval(int index, int i)
	{
		return comp.getmultival(index, i);
	}
	void setmulti(vector<vector<int>>& edges, int grainsize, int edgefact, int starts)
	{
		comp.setmulti(edges, grainsize, edgefact, starts);
	}
	void multival(vector<int> starts)
	{
		comp.multival(starts, false);
	}

	void bfs_search(vector<int> starts)
	{
		comp.bfs_search_act(starts);
	}
	int getmultival(int index, int i)
	{

		return comp.getmultival(index, i);
	}
	int getnum()
	{

		return comp.getnum();
	}

};

void parallel_edge_gen(vector<int>::iterator pedges, vector<vector<int>>* nodes, int size);
#endif