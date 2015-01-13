
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
	vector<std::shared_future<void>> threads;
	
	~ThreadBag()
	{
	}
	void add(std::vector<pair<Edge, Edge>> holder, int loc, vector<SubGraph*> neighbors)
	{
		for (int i = 0; i < holder.size(); ++i)
		{
			pair<Edge, Edge> temp = holder[i];
			Edge first = temp.first;
			Edge second = temp.second;
			threads.push_back(
				std::shared_future<void>(std::async(
				&bfs_search_act, first.first, first.second, loc, second.first, neighbors[second.second])
				)
				);
		}
	}
	void collect()
	{
		int count = 0;
		for (int j = 0; j < threads.size(); ++j)
		{
			if (threads[j].valid())
				threads[j].get();
			++count;
		}
		collected = count;
		//threads = vector<std::shared_future<void>>();
	}
	int collected = 0;
};

struct EdgeBag
{
	vector<pair<Edge, int>> v;
	int spot = 0;
	static void locker(bool lock)
	{
		static std::mutex m;
		if (lock)
			while (!m.try_lock());
		else
			m.unlock();
			
	}

	void add(pair<Edge,int> t)
	{
		locker(true);
		v.push_back(t);
		locker(false);
	}
	pair<Edge,int> get()
	{
		locker(true);
		pair<Edge,int> index;
		index = v[spot];
		++spot;
		locker(false);
		return index;
	}
	int size()
	{
		locker(true);
		int retval = 0;
		retval = v.size()-spot;
		locker(false);
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
	ThreadBag b[64];

	bool lock(int loc, bool unlock)
	{
		static std::mutex m[64];
		if (!unlock)
		{
			if (!m[loc].try_lock())
				return false;
			return true;
		}
		else
			m[loc].unlock();
		return true;
	}
	EdgeBag v[64];
	void search(int loc)
	{
		
		property_map < BoolGraph, vertex_index_t >::type
			index_map = get(vertex_index, g);
		int parent;
		vector<pair<Edge, Edge>> holder;
		holder.reserve(128);
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
					holder.push_back(pair<Edge, Edge>(Edge(ind, parent), Edge(val.second + 1, spot)));
					if (holder.size() >= grainsize)
					{
						b[loc].add(holder, loc, neighbors);
						holder = vector<pair<Edge, Edge>>();
					}
					continue;
				}
				v[loc].add(pair<Edge, int>(Edge(ind, parent), val.second + 1));
			}
		}
		lock(loc, true);
		if (holder.size() > 0)
		b[loc].add(holder, loc, neighbors);
	}
	bool bfs_init(int loc)
	{
		if (v[loc].size() <= 0)
			return true;
		if (!lock(loc, false))
		{
			return false;
		}
		search(loc);
		return false;
	}
	void bfs_search(int index, int parent, int loc, int dist)
	{
		if (!lock(loc, false))
		{
			if (name[parent][loc].second > dist - 1)
			{
				name[parent][loc].second = dist-1;
			}
			if (name[index][loc].second <= dist)
			{
				return;
			}
			v[loc].add(pair<Edge,int>(Edge(index, parent), dist));
			return;
		}
		if (name[index][loc].second <= dist)
		{
			if (name[parent][loc].second > dist - 1 && dist != 0)
			{
				name[parent][loc].second = dist - 1;
			}

			lock(loc, true);
			return;
		}
		if (name[parent][loc].second > dist - 1 && dist != 0)
		{
			name[parent][loc].second = dist - 1;
		}

		v[loc].add(pair<Edge, int>(Edge(index, parent), dist));
		search(loc);
	};

	int getnum()
	{
		return num_vertices(g);
	}
	int getval(int i, int index)
	{
		return name[i][index].first;
	}
	bool nothreads(int loc)
	{
		if (b[loc].threads.size()-b[loc].collected <= 0)
		{
			if (lock(loc, false))
			{
				lock(loc, true);
				return bfs_init(loc);
			}
			return false;
		}
		if (lock(loc, false))
		{
			lock(loc, true);
			b[loc].collect();
			return false;
		}
		return false;
	}
	property_map<BoolGraph, bool_name_t>::type
		name;
	BoolGraph g;
	int grainsize = 9999;
	int edgefact = 16;
	vector<int> parts;
	vector<SubGraph*> neighbors;
	int part;
	bool done = false;
};

static void bfs_search_act(int index, int parent, int loc, int dist, SubGraph* g)
{
	g->bfs_search(index, parent, loc, dist);
}
struct MultiComponent
{
	void set(vector<int> parts, vector<vector<int>> nodes, int size, int edge, int starts, int nparts)
	{
		graphs = vector<SubGraph>(nparts);
		grainsize = size;
		edgefact = edge;
		vertexmap = parts;
		for (int i = 0; i < nparts; ++i)
		{
			graphs[i].set(nodes, size, edge, starts, i, vertexmap);
		}
		neighbors = vector<SubGraph*>(nparts);
		for (int i = 0; i < graphs.size(); ++i)
		{
			neighbors[i] = &graphs[i];
		}
		for (int i = 0; i < graphs.size(); ++i)
		{
			graphs[i].neighbors = neighbors;
		}
		parents = vector<vector<int>>(graphs.front().getnum(), vector<int>(starts, -1));
	};
	static void run(int start, int size, int i, vector<SubGraph>* graphs, vector<vector<int>>* parents, vector<int>* vertexmap)
	{
		int index = (*vertexmap)[start];
		(*graphs)[index].bfs_search(start, start, i, 0);
		bool waiter = true;
		while (waiter)
		{
			waiter = false;
			for (int j = 0; j < graphs->size(); ++j)
			{
				if (!(*graphs)[j].nothreads(i))
				{
					waiter = true;
					break;
				}
			}
		}
		for (int j = 0; j < vertexmap->size(); ++j)
		{
			index = (*vertexmap)[j];
			(*parents)[j][i] = (*graphs)[index].getval(j, i);
		}
	}
	void search(vector<int> starts)
	{
		int size = grainsize;
		vector<std::thread> runs(starts.size());
		for (int i = 0; i < starts.size(); ++i)
		{
			int start = starts[i];
			runs[i] = std::thread(&run, start, size, i, &graphs, &parents, &vertexmap);
		}
		for (int i = 0; i < starts.size(); ++i)
			runs[i].join();
	}
	int getmultival(int index, int i)
	{
		return parents[index][i];
	}
	vector<SubGraph> graphs;
	vector<int> vertexmap;
	vector<vector<int>> parents;
	int grainsize = 9999;
	int edgefact = 16;
	vector<SubGraph*> neighbors;
};

