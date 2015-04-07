
#include "headers.hpp"

struct bool_name_t {
	typedef boost::vertex_property_tag kind;
};

typedef boost::property<bool_name_t, std::vector<Edge> > BoolColor; //parent, color, partition, distance
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	BoolColor> BoolGraph;

typedef struct SubGraph;
void bfs_search_act(int index, int parent, int loc, int dist, SubGraph* g);

typedef hpx::lcos::local::spinlock locker;

struct ThreadBag
{
	vector<hpx::shared_future<void>> threads;
	vector<SubGraph*> neighbors;
	vector<pair<Edge, Edge>> edges;
	int loc = -1;
	int count = 0;
	int part = -1;
	int grainsize;
	ThreadBag(int gs)
	{
		grainsize = gs;
		edges.reserve(grainsize);
	}
	~ThreadBag()
	{
		collect();
	}
	void collect()
	{

		for (int i = 0; i < edges.size(); ++i)
		{
			pair<Edge, Edge> temp = edges[i];
			Edge first = temp.first;
			Edge second = temp.second;
			hpx::async(
				&bfs_search_act, first.first, first.second, loc, second.first, neighbors[second.second]);
		}
	}
	void add(pair<Edge, Edge> holder)
	{
		edges.push_back(holder);
		if (edges.size() >= grainsize)
		{
			for (int i = 0; i < edges.size(); ++i)
			{
				pair<Edge, Edge> temp = edges[i];
				Edge first = temp.first;
				Edge second = temp.second;
				hpx::async(
					&bfs_search_act, first.first, first.second, loc, second.first, neighbors[second.second]);
			}
			edges = vector<pair<Edge, Edge>>();
			edges.reserve(grainsize);
		}
	}
};

struct EdgeBag
{
	mutable locker l;
	vector<pair<Edge, int>> v;
	int spot = 0;

	void add(pair<Edge,int> t)
	{
		locker::scoped_lock m(l);
		v.push_back(t);
	}
	pair<Edge,int> get()
	{
		locker::scoped_lock m(l);
		pair<Edge,int> index;
		index = v[spot];
		++spot;
		return index;
	}
	int size()
	{
		locker::scoped_lock m(l);
		int retval = 0;
		retval = v.size() - spot;
		return retval;
		
	}
};

struct SubGraph
{

	void multireset(int starts)
	{
		for (int i = 0; i < num_vertices(g); ++i)
		{
			name[i] = vector<Edge>(starts, Edge(-1, 99999999));
		}
	}

	void set(vector<vector<int>> nodes, int size, int edge, int starts, int index, vector<int> map)
	{
		part = index;
		parts = map;
		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int j = 0; j < nodes[i].size(); ++j)
			{
				int val = nodes[i][j];
				if (parts[i] == part || parts[val] == part)
					add_edge(i, val, g);
			}
		}
		grainsize = size;
		edgefact = edge;
		
		name = get(bool_name_t(), g);
		multireset(starts);
	}

	EdgeBag v[64];
	void search(pair<Edge, int> sample, const int loc)
	{
		v[loc].add(sample);
		if (!l[loc].try_lock())
			return;
		ThreadBag b(grainsize);
		b.neighbors = neighbors;
		b.loc = loc;
		property_map < BoolGraph, vertex_index_t >::type
			index_map = get(vertex_index, g);
		int parent;
		while (v[loc].size() > 0)
		{
			pair<Edge, int> val = v[loc].get();
			parent = val.first.first;
			if (val.second > name[parent][loc].second)
				continue;
			name[parent][loc] = Edge(val.first.second, val.second);
			graph_traits < BoolGraph >::adjacency_iterator ai, a_end;
			for (boost::tie(ai, a_end) = adjacent_vertices(parent, g); ai != a_end; ++ai)
			{
				int ind = get(index_map, *ai);
				if (name[ind][loc].second <= val.second + 1)
					continue;

				name[ind][loc] = Edge(parent, val.second + 1);
				int spot = parts[ind];
				if (spot != part)
				{
					b.add(pair<Edge, Edge>(Edge(ind, parent), Edge(val.second + 1, spot)));
					continue;
				}
				v[loc].add(pair<Edge, int>(Edge(ind, parent), val.second + 1));
			}
		}
		l[loc].unlock();
	}
	void bfs_search(int index, int parent, int loc, int dist)
	{
		if (name[index][loc].second <= dist)
		{
			if (name[parent][loc].second > dist - 1)
			{
				name[parent][loc].second = dist - 1;
			}
			return;
		}
		if (name[parent][loc].second > dist - 1 && dist != 0)
		{
			name[parent][loc].second = dist - 1;
		}

		search(pair<Edge, int>(Edge(index, parent), dist), loc);
	};

	int getnum()
	{
		return num_vertices(g);
	}
	int getval(int i, int index)
	{
		return name[i][index].first;
	}
	property_map<BoolGraph, bool_name_t>::type
		name;
	BoolGraph g;
	int grainsize = 9999;
	int edgefact = 16;
	vector<int> parts;
	vector<SubGraph*> neighbors;
	mutable locker l[64];
	int part;
};
void bfs_search_act(int index, int parent, int loc, int dist, SubGraph* g)
{
	g->bfs_search(index, parent, loc, dist);
}

struct MultiComponent
{
	void set(vector<int> parts, vector<vector<int>> nodes, int size, int edge, int starts, int nparts)
	{
		//graphs = vector<SubGraph>(nparts, SubGraph());
		graphs.resize(nparts);
		grainsize = size;
		edgefact = edge;
		vertexmap = parts;
		neighbors = vector<SubGraph*>(nparts);
		for (int i = 0; i < nparts; ++i)
		{
			graphs[i].set(nodes, size, edge, starts, i, vertexmap);
			neighbors[i] = &graphs[i];
		}
		for (auto it = graphs.begin(); it != graphs.end(); ++it)
		{
			it->neighbors = neighbors;
		}
	};
	static void run(int start, int size, int i, vector<SubGraph*> graphs, vector<int>* vertexmap)
	{
		int index = (*vertexmap)[start];
		(graphs)[index]->bfs_search(start, start, i, 0);
	}
	void search(vector<int> starts)
	{
		int size = grainsize;
		vector<hpx::future<void>> runs(starts.size());
		for (int i = 0; i < starts.size(); ++i)
		{
			int start = starts[i];
			runs[i] = hpx::async(&run, start, size, i, neighbors, &vertexmap);
		}
		hpx::wait_all(runs);
		hpx::finalize();
	}
	int getmultival(int index, int i)
	{
		int start = vertexmap[index];
		return graphs[start].getval(index,i);
	}
	deque<SubGraph> graphs;
	vector<int> vertexmap;
	//vector<vector<int>> parents;
	int grainsize = 9999;
	int edgefact = 16;
	vector<SubGraph*> neighbors;
};

