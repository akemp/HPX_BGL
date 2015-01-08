//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
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
				if (parts[i] == part || parts[val] == part)
					add_edge(i, val, g);
			}
		}
		grainsize = size;
		edgefact = edge;

		name = get(bool_name_t(), g);
		multireset(starts);
	}
	vector<pair<int,int>> bfs_search(vector<pair<int,int>> indices, int loc)
	{
		vector<pair<int,int>> q;
		property_map < BoolGraph, vertex_index_t >::type
			index_map = get(vertex_index, g);
		for (int i = 0; i < indices.size(); ++i)
		{
			int parent = (indices)[i].first;
			name[parent][loc] = true;
			if (parts[parent] != part)
				continue;
			graph_traits < BoolGraph >::adjacency_iterator ai, a_end;

			for (boost::tie(ai, a_end) = adjacent_vertices(parent, g); ai != a_end; ++ai)
			{
				int ind = get(index_map, *ai);
				if (name[ind][loc])
					continue;
				name[ind][loc] = true;
				q.push_back(pair<int, int>(ind, parent));
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
		parents = vector<vector<int>>(graphs.front().getnum(), vector<int>(starts, -1));
	};
	static vector<pair<int, int>> bfs_search_act(SubGraph* graph, vector<pair<int, int>> starter, int i)
	{
		return graph->bfs_search(starter, i);
	}
	void search(vector<int> starts)
	{
		int size = grainsize;
		for (int i = 0; i < starts.size(); ++i)
		{
			int start = starts[i];
			vector<pair<int,int>> starter;
			starter.push_back(pair<int,int>(start,start));
			while (starter.size() > 0)
			{
				vector<hpx::future<vector<pair<int, int>>>> futs;

				int spot = 0;
				for (vector<pair<int, int>>::iterator it = starter.begin(); it < starter.end(); it += size)
				{
					int last = spot;
					spot += size;
					vector<pair<int, int>> input;
					if (spot < starter.size())
						input = vector<pair<int, int>>(it, it + size);
					else
						input = vector<pair<int, int>>(it, it + starter.size() - last);
					for (int j = 0; j < graphs.size(); ++j)
					{
						futs.push_back(hpx::async(&bfs_search_act, &graphs[j], input, i));
					}
					if (spot >= starter.size())
					break;
				}
				starter.clear();
				for (int j = 0; j < futs.size(); ++j)
				{
					vector<pair<int, int>> child = futs[j].get();;
					for (int k = 0; k < child.size(); ++k)
					{
						pair<int, int> sample = child[k];
						if (parents[sample.first][i] < 0)
						{
							parents[sample.first][i] = sample.second;
							starter.push_back(sample);
						}
					}
				}
			}
		}
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
};


///////////////////////////////////////////////////////////////////////////////
int main()
{

	using namespace std;

    uint64_t nnodes;
    uint64_t scale;
	cout << "Enter scale (nodes=2^input): ";
	cin >> scale;
    uint64_t edgefactor;
	cout << "Enter edgefactor: ";
    cin >> edgefactor;
	int grainsize;
#ifdef CUSTOMGRAIN
	cout << "Enter grainsize: ";
	cin >> grainsize;
#else
	grainsize = 128;
#endif

	int searches;
	cout << "Enter searches: ";
	cin >> searches;
	int acctest = 1;
	//cout << "Run accuracy tests (0 for no, 1 for yes)?";
	//cin >> acctest;
    nnodes = (uint64_t)(1) << scale;
	//
	int nparts = 10;
	cout << "Enter partitions: ";
	cin >> nparts;


	//Edges edges;
	vector<vector<int>> nodes(
		nnodes
		);
    uint64_t ind = (uint64_t)(edgefactor) << scale;
	vector<int> pedges(ind);
	boost::random::mt19937 rng;
	boost::random::uniform_int_distribution<> randnodes(0, nnodes-1);
	for (int i = 0; i < pedges.size(); ++i)
		pedges[i] = randnodes(rng);
	cout << "Random range generated. Making edgelist.\n";
	{
		parallel_edge_gen(pedges.begin(),&nodes, ind);
	}

	cout << "Edgelist generated. Creating partitions.\n";
	vector<int> xadj, adjncy;
	toCSR(nodes, xadj, adjncy);
	vector<int> parts(nodes.size());
	get_parts(xadj, adjncy, parts, nparts);
	cout << "Running tests.\n";


	vector<int> starts(searches);


	vector<vector<pair<int, int>>> counts;
	cout << "Setting up serial values for testing.\n";
	counts = vector<vector<pair<int, int>>>(starts.size(), vector<pair<int, int>>(64, pair<int, int>(-1, 0)));
	cout << "Data structures set. Running serial search.\n";
	int nverts;
	{
		graph_manager hw;// = graph_manager::create(hpx::find_here());
		hw.setmulti(nodes, grainsize, ind, starts.size());
		nverts = hw.getnum();
		randnodes = boost::random::uniform_int_distribution<>(0, nverts - 1);
		for (int j = 0; j < starts.size(); ++j)
		{
			starts[j] = randnodes(rng);
		}

        hpx::util::high_resolution_timer t1;
        hw.bfs_search(starts);
        double elapsed1 = t1.elapsed();
        cout << elapsed1 << "s search time for serial\n";
		for (int j = 0; j < starts.size(); ++j)
		{
			for (int i = 0; i < counts[j].size(); ++i)
			{
				int sample = randnodes(rng);
				counts[j][i].first = sample;
				int count = 0;
				while (sample != starts[j])
				{
					count++;
					if (sample == -1)
					{
						count = -1;
                        cout << "Degenerate vertex in sample " << j << "-" << i << endl;
						break;
					}
					//cout << sample << "-" << sub.pennants[sample].dist << " ";
					sample = hw.getmultival(sample, j);
				}
				counts[j][i].second = count;
				//cout << endl;
			}
		}
	}
	if (acctest != 0)
	{
		cout << "Running accuracy tests.\n";
		

			hpx::util::high_resolution_timer t;
			graph_manager hw;// = graph_manager::create(hpx::find_here());
		  hw.set(nodes, grainsize, ind, starts.size());
		hpx::util::high_resolution_timer t1;
		hw.pbfs_search(starts);
		double elapsed = t.elapsed();
		double elapsed1 = t1.elapsed();
		cout << elapsed << "s for parallel component\n";
		cout << elapsed1 << "s search time for parallel component\n";
		for (int j = 0; j < starts.size(); ++j)
		{
		  
			for (int i = 0; i < counts[j].size(); ++i)
			{
				int sample = counts[j][i].first;
				int count = 0;
				while (sample != starts[j])
				{
					count++;
					if (sample == -1)
					{
						count = -1;
						break;
					}
					//cout << sample << "-" << sub.pennants[sample].dist << " ";
					sample = hw.getval(sample, j);
					if (counts[j][i].second < count)
						break;
				}
				if (counts[j][i].second != count)
					cout << "Counts not equal! " << count << " for bfs != " << counts[j][i].second << " for pbfs!\n";
				//cout << endl;

		  }
	  }
	  {
		  hpx::util::high_resolution_timer t;
		  t.restart();
		  graph_manager hw;// = graph_manager::create(hpx::find_here());
		  hw.setmulti(nodes, grainsize, ind, starts.size());
		  hpx::util::high_resolution_timer t1;
		  hw.multival(starts);

		  double elapsed1 = t1.elapsed();
		  double elapsed = t.elapsed();
		  cout << elapsed << "s for highly parallel\n";
		  cout << elapsed1 << "s search time for highly parallel\n";


		  for (int j = 0; j < starts.size(); ++j)
		  {
			  for (int i = 0; i < counts[j].size(); ++i)
			  {
				  int sample = counts[j][i].first;
				  int count = 0;
				  while (sample != starts[j])
				  {
					  count++;
					  if (sample == -1)
					  {
						  count = -1;
						  break;
					  }
					  //cout << sample << "-" << sub.pennants[sample].dist << " ";
					  sample = hw.getmultival(sample, j);
					  if (counts[j][i].second < count)
						  break;
				  }
				  if (counts[j][i].second != count)
					  cout << "Counts not equal! " << count << " for bfs != " << counts[j][i].second << " for pbfs!\n";
				  //cout << endl;
			  }
		  }
	  }

		{
			hpx::util::high_resolution_timer t;
			t.restart();
			//MultiComponent(vector<int> parts, vector<vector<int>> nodes, int size, int edge, int starts, int nparts)
			MultiComponent hw;// = graph_manager::create(hpx::find_here());

			//MultiComponent(vector<int> parts, vector<vector<int>> nodes, int size, int edge, int starts, int nparts)
			hw.set(parts, nodes, grainsize, ind, starts.size(), nparts);

			hpx::util::high_resolution_timer t1;
			{
				hw.search(starts);
			}
			double elapsed = t.elapsed();
			double elapsed1 = t1.elapsed();
			cout << elapsed << "s for partitioned\n";
			cout << elapsed1 << "s search time for partitioned\n";
			for (int j = 0; j < starts.size(); ++j)
			{
				for (int i = 0; i < counts[j].size(); ++i)
				{
					int sample = counts[j][i].first;
					int count = 0;
					while (sample != starts[j])
					{
						count++;
						if (sample == -1)
						{
							count = -1;
							break;
						}
						//cout << sample << "-" << sub.pennants[sample].dist << " ";
						sample = hw.getmultival(sample, j);
						if (counts[j][i].second < count)
							break;
					}
					if (counts[j][i].second != count)
						cout << "Counts not equal! " << count << " for pbfs != " << counts[j][i].second << " for bfs!\n";
					//cout << endl;
				}
			}
		}
		cout << "Accuracy tests complete.\n";
	}
	cout << "Serial comparison:\n";
	{
		hpx::util::high_resolution_timer t;
		t.restart();
		graph_manager hw;// = graph_manager::create(hpx::find_here());
		hw.setmulti(nodes, grainsize, ind, starts.size());

		hpx::util::high_resolution_timer t1;
		{
			hw.bfs_search(starts);
		}
		double elapsed = t.elapsed();
		double elapsed1 = t1.elapsed();
		cout << elapsed << "s for serial\n";
		cout << elapsed1 << "s search time for serial\n";
	}

	int s;
	cin >> s;
	return 0;
}