
#include "headers.hpp"

struct bool_name_t {
	typedef boost::vertex_property_tag kind;
};

typedef boost::property<bool_name_t, std::vector<int> > BoolColor; //parent, color, partition, distance
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	BoolColor> BoolGraph;

struct ThreadBag
{
	vector<std::shared_future<void>> threads;

	~ThreadBag()
	{
		for (int j = 0; j < threads.size(); ++j)
		{
			threads[j].wait();
		}
	}
	void add(std::shared_future<void> t)
	{
		static std::mutex m;
		while (!m.try_lock());
		threads.push_back(t);
		m.unlock();
	}
	void collect()
	{
		for (int j = 0; j < threads.size(); ++j)
		{
			if (threads[j].valid())
				threads[j].wait();
		}
		threads = vector<std::shared_future<void>>();
	}
};

struct EdgeBag
{
	vector<pair<Edge, int>> v;

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
		pair<Edge,int> index;
		locker(true);
		if (v.size() <= 0)
		{
			locker(false);
			return pair<Edge, int>(Edge(-1, -1), -1);
		}
		index = v.front();
		v.erase(v.begin());
		locker(false);
		return index;
	}
	int size()
	{
		int retval = 0;
		locker(true);
		retval = v.size();
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
			name[i] = vector<int>(starts, -1);
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
		dists = vector<vector<int>>(nodes.size(), vector<int>(starts, 99999999999999));
		multireset(starts);
	}
	static void bfs_search_act(int index, int parent, int loc, int dist, SubGraph* g)
	{
		g->bfs_search(index, parent, loc, dist);
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
	bool bfs_init(int loc)
	{
		if (v[loc].size() <= 0)
			return true;
		pair<Edge, int> val = v[loc].get();
		if (val.second < 0)
			return true;
		bfs_search(val.first.first, val.first.second, loc, val.second);
		return false;
	}
	void bfs_search(int index, int parent, int loc, int dist)
	{
		if (!lock(loc, false))
		{
			//std::this_thread::sleep_for(std::chrono::microseconds(1000));
			if (dists[index][loc] <= dist)
			{
				return;
			}

			v[loc].add(pair<Edge,int>(Edge(index, parent), dist));
			return;
		}
		if (dists[index][loc] <= dist)
		{
			lock(loc, true);
			return;
		}
		property_map < BoolGraph, vertex_index_t >::type
		index_map = get(vertex_index, g);
		v[loc].add(pair<Edge, int>(Edge(index, parent), dist));
		while (v[loc].size() > 0)
		{
			pair<Edge, int> val = v[loc].get();
			Edge sample = val.first;
			parent = sample.first;
			if (val.second > dists[parent][loc])
				continue;
			dists[parent][loc] = val.second;
			name[parent][loc] = sample.second;
			graph_traits < BoolGraph >::adjacency_iterator ai, a_end;
			for (boost::tie(ai, a_end) = adjacent_vertices(parent, g); ai != a_end; ++ai)
			{
				int ind = get(index_map, *ai);
				if (dists[ind][loc] <= dists[parent][loc] + 1)
					continue;

				dists[ind][loc] = dists[parent][loc] + 1;
				name[ind][loc] = parent;
				int spot = parts[ind];
				if (spot != part)
				{
					sample = Edge(ind, parent);
					int indexer = parts[sample.first];
					b[loc].add(
						std::shared_future<void>(std::async(
						&bfs_search_act, sample.first, sample.second, loc, dists[sample.first][loc], neighbors[indexer])
						)
						);
					continue;
				}
				v[loc].add(pair<Edge, int>(Edge(ind, parent), dists[parent][loc] + 1));
			}
		}
		lock(loc, true);
	};

	int getnum()
	{
		return num_vertices(g);
	}
	int getval(int i, int index)
	{
		return name[i][index];
	}
	bool nothreads(int loc)
	{
		if (b[loc].threads.size() <= 0)
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
	vector<vector<int>> dists;
};

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

