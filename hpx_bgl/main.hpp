
#include "headers.hpp"

struct bool_name_t {
	typedef boost::vertex_property_tag kind;
};

typedef boost::property<bool_name_t, std::vector<bool> > BoolColor; //parent, color, partition, distance
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	BoolColor> BoolGraph;

struct SubGraph
{
	~SubGraph()
	{

	}

	void multireset(int starts)
	{
		for (int i = 0; i < num_vertices(g); ++i)
		{
			name[i] = vector<bool>(starts, false);
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
				//if (parts[i] == part || parts[val] == part)
					add_edge(i, val, g);
			}
		}
		grainsize = size;
		edgefact = edge;

		name = get(bool_name_t(), g);
		multireset(starts);
	}
	vector<pair<int, int>> bfs_search(int index, int parent, int loc)
	{
		vector<pair<int, int>> q;
		if (name[index][loc])
			return q;
		property_map < BoolGraph, vertex_index_t >::type
			index_map = get(vertex_index, g);
		q.push_back(pair<int,int>(index,parent));
		{
			for (int i = 0; i < q.size(); ++i)
			{
				parent = q[i].first;
				name[parent][loc] = true;
				graph_traits < BoolGraph >::adjacency_iterator ai, a_end;

				for (boost::tie(ai, a_end) = adjacent_vertices(parent, g); ai != a_end; ++ai)
				{
					int ind = get(index_map, *ai);
					if (name[ind][loc])
						continue;
					int spot = parts[ind];
					if (spot != part)
					{
						vector<pair<int, int>> temp = neighbors[spot]->bfs_search(ind, parent, loc);
						q.insert(q.end(), temp.begin(), temp.end());
						continue;
					}
					q.push_back(pair<int, int>(ind, parent));
				}
			}
		}
		return q;
	};

	int getnum()
	{
		return num_vertices(g);
	}
	property_map<BoolGraph, bool_name_t>::type
		name;
	BoolGraph g;
	int grainsize = 9999;
	int edgefact = 16;
	vector<int> parts;
	vector<SubGraph*> neighbors;
	int part;
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
	static void run(int start, int size, int i, vector<SubGraph>* graphs, vector<vector<int>>* parents)
	{

		vector<pair<int, int>> child = (graphs)->front().bfs_search(start, start, i);
		for (int k = 0; k < child.size(); ++k)
		{
			pair<int, int> sample = child[k];
			if ((*parents)[sample.first][i] < 0)
			{
				(*parents)[sample.first][i] = sample.second;
			}
		}
		/*
		while (starter.size() > 0)
		{
			vector<hpx::future<vector<pair<int, int>>>> futs;

			int spot = 0;
			vector < vector<pair<int, int>>> inputs;
			int counter = 0;
			
			for (vector<pair<int, int>>::iterator it = starter.begin(); it < starter.end(); it += size)
			{
				int last = spot;
				spot += size;

				if (spot < starter.size())
					inputs.push_back(vector<pair<int, int>>(it, it + size));
				else
					inputs.push_back(vector<pair<int, int>>(it, it + starter.size() - last));
				for (int j = 0; j < graphs->size(); ++j)
				{
					futs.push_back(hpx::async(&bfs_search_act, &(*graphs)[j], inputs[counter], i));
				}
				if (spot >= starter.size())
					break;
				++counter;
			}
			starter.clear();
			for (int j = 0; j < futs.size(); ++j)
			{
				vector<pair<int, int>> child = futs[j].get();;
				for (int k = 0; k < child.size(); ++k)
				{
					pair<int, int> sample = child[k];
					if ((*parents)[sample.first][i] < 0)
					{
						(*parents)[sample.first][i] = sample.second;
						starter.push_back(sample);
					}
				}
			}
		}*/
	}
	void search(vector<int> starts)
	{
		int size = grainsize;
		vector<hpx::thread> runs(starts.size());
		for (int i = 0; i < starts.size(); ++i)
		{
			int start = starts[i];
			runs[i] = hpx::thread(&run, start, size, i, &graphs, &parents);
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

