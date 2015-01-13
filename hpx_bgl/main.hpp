
#include "headers.hpp"

struct bool_name_t {
	typedef boost::vertex_property_tag kind;
};

typedef boost::property<bool_name_t, std::vector<Edge> > BoolColor; //parent, color, partition, distance
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	BoolColor> BoolGraph;

typedef struct SubGraph;
static void bfs_search_act(int index, int parent, int loc, int dist, SubGraph* g);

struct ThreadBag
{
	vector<hpx::shared_future<void>> threads;
	vector<SubGraph*> neighbors;
	~ThreadBag()
	{
		collect();
	}
	void add(pair<Edge, Edge> holder, int loc)
	{
		//m.lock();
		//for (int i = 0; i < holder.size(); ++i)
		{
			pair<Edge, Edge> temp = holder;
			Edge first = temp.first;
			Edge second = temp.second;
			threads.push_back(
				hpx::shared_future<void>(hpx::async(
				&bfs_search_act, first.first, first.second, loc, second.first, neighbors[second.second])
				)
				);
		}
		//m.unlock();
	}
	void collect()
	{
		hpx::wait_all(threads);
	}
};

struct EdgeBag
{
	vector<pair<Edge, int>> v;
	int spot = 0;

	void add(pair<Edge,int> t)
	{
		v.push_back(t);
	}
	pair<Edge,int> get()
	{
		pair<Edge,int> index;
		index = v[spot];
		++spot;
		return index;
	}
	int size()
	{
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

	bool lock(int loc, bool unlock)
	{
		if (!unlock)
		{
			m[loc].lock();
		}
		else
			m[loc].unlock();
		return true;
	}
	void search(pair<Edge, int> sample, int loc)
	{
		lock(loc, false);
		ThreadBag b;
		b.neighbors = neighbors;
		EdgeBag v;
		v.add(sample);
		property_map < BoolGraph, vertex_index_t >::type
			index_map = get(vertex_index, g);
		int parent;
		while (v.size() > 0)
		{
			pair<Edge, int> val = v.get();
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
					//holder.push_back(pair<Edge, Edge>(Edge(ind, parent), Edge(val.second + 1, spot)));
					//if (holder.size() >= grainsize)
					{
						b.add(pair<Edge, Edge>(Edge(ind, parent), Edge(val.second + 1, spot)), loc);
						//holder = vector<pair<Edge, Edge>>();
					}
					continue;
				}
				v.add(pair<Edge, int>(Edge(ind, parent), val.second + 1));
			}
		}
		lock(loc, true);
	}
	void bfs_search(int index, int parent, int loc, int dist)
	{
		if (name[index][loc].second <= dist)
		{
			if (name[parent][loc].second > dist - 1 && dist != 0)
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
	int part;
	hpx::lcos::local::mutex m[64];
};

static void bfs_search_act(int index, int parent, int loc, int dist, SubGraph* g)
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
		parents = vector<vector<int>>(graphs.front().getnum(), vector<int>(starts, -1));
	};
	static void run(int start, int size, int i, vector<SubGraph*> graphs, vector<vector<int>>* parents, vector<int>* vertexmap)
	{
		int index = (*vertexmap)[start];
		(graphs)[index]->bfs_search(start, start, i, 0);
		for (int j = 0; j < vertexmap->size(); ++j)
		{
			index = (*vertexmap)[j];
			(*parents)[j][i] = (graphs)[index]->getval(j, i);
		}
	}
	void search(vector<int> starts)
	{
		int size = grainsize;
		vector<hpx::future<void>> runs(starts.size());
		for (int i = 0; i < starts.size(); ++i)
		{
			int start = starts[i];
			runs[i] = hpx::async(&run, start, size, i, neighbors, &parents, &vertexmap);
		}
		hpx::wait_all(runs);
	}
	int getmultival(int index, int i)
	{
		return parents[index][i];
	}
	deque<SubGraph> graphs;
	vector<int> vertexmap;
	vector<vector<int>> parents;
	int grainsize = 9999;
	int edgefact = 16;
	vector<SubGraph*> neighbors;
};

