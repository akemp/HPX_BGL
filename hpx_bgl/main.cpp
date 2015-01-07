//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include "headers.hpp"

struct SubGraph
{
	~SubGraph()
	{
	}
	static void parreset(PairGraph* g, int start, int size, int toggled)
	{
		property_map<PairGraph, pair_name_t>::type
			name = get(pair_name_t(), *g);
		for (int i = start; i < start + size; ++i)
		{
			name[i] = vector<pair<int,int>>(toggled, pair<int,int>(-1,9999999));
		}
	}

	void reset(int toggled)
	{
		int adder = grainsize;

		vector<std::thread> futs;
		PairGraph* ptr = &g;

		for (int i = 0; i < num_vertices(g); i += adder)
		{
			if (i + adder < num_vertices(g))
			{
				futs.push_back(std::thread(&parreset, ptr, i, adder, toggled));
			}
			else
			{
				futs.push_back(std::thread(&parreset, ptr, i, num_vertices(g) - i, toggled));
			}
		}
		for (int i = 0; i < futs.size(); ++i)
			futs[i].join();

	}

	int getval(int i, int index)
	{
		return name[i][index].first;
	}
	void multireset(int starts)
	{
		for (int i = 0; i < num_vertices(g); ++i)
		{
			name[i] = vector<pair<int,int>>(starts, pair<int,int>(-1,9999999));
		}
	}

	void set(vector<vector<int>> nodes, int size, int edge, int starts, int index, vector<int> map, vector<SubGraph*> n)
	{
		part = index;
		parts = map;
		grainsize = size;
		edgefact = edge;
		neighbors = n;
		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int j = 0; j < nodes[i].size(); ++j)
			{
				int val = nodes[i][j];
				if (parts[i] == part || parts[val] == part)
					add_edge(i, val, g);
			}
		}
		name = get(pair_name_t(), g);
		multireset(starts);
	}

	int getmultival(int index, int i)
	{
		return name[index][i].first;
	}

	static void bfs_partition(int start, int loc, SubGraph* g, int parent, int dist)
	{
		g->bfs_search(start, loc, parent, dist);
	}
	void bfs_search(int index, int loc, int parent, int dist)
	{
		property_map < PairGraph, vertex_index_t >::type
			index_map = get(vertex_index, g);
		//pennants = std::vector <int>(num_vertices(g), -1);
		if (dist >= name[index][loc].second)
		{
			return;
		}

		name[parent][loc].second = dist - 1;
		
		name[index][loc].first = parent;
		name[index][loc].second = dist;

		vector<int> q;
		q.push_back(index);
		int spot = 0;
		while (spot < q.size())
		{
			int ind = q[spot];
			++spot;
			parent = ind;
			graph_traits < PairGraph >::adjacency_iterator ai, a_end;

			for (boost::tie(ai, a_end) = adjacent_vertices(ind, g); ai != a_end; ++ai)
			{
				ind = get(index_map, *ai);
				int dist2 = name[parent][loc].second + 1;
				if (dist2 < name[ind][loc].second)
				{
					name[ind][loc].second = dist2;
					name[ind][loc].first = parent;
					int sampart = parts[ind];
					if (sampart != part)
					{
						bfs_partition(ind, loc, neighbors[sampart], parent, dist2);
						continue;
					}
					q.push_back(ind);
				}
			}

		}
	};

	int getnum()
	{
		return num_vertices(g);
	}

	property_map<PairGraph, pair_name_t>::type
		name;
	PairGraph g;
	int grainsize = 9999;
	int edgefact = 16;
	vector<SubGraph*> neighbors;
	vector<int> parts;
	int part;
	vector<shared_future<void>> threads;
	int numthreads = 0;
};

struct MultiComponent
{
	void set(vector<int> parts, vector<vector<int>> nodes, int size, int edge, int starts, int nparts)
	{
		graphs = vector<SubGraph>(nparts);
		neighbors = vector<SubGraph*>(nparts);
		grainsize = size;
		edgefact = edge;
		vertexmap = parts;
		//void set(vector<vector<int>> nodes, int size, int edge,
		//int starts, int partno, vector<int> map, int index, vector<SubGraph*> n)
		for (int i = 0; i < graphs.size(); ++i)
		{
			neighbors[i] = (&(graphs[i]));
		}
		for (int i = 0; i < nparts; ++i)
		{
			graphs[i].set(nodes, grainsize, edgefact, starts, i, vertexmap, neighbors);
		}
	};
	static void search_act(int start, int i, SubGraph* g)
	{
		g->bfs_search(start, i, start, 0);
	}
	void search(vector<int> starts)
	{
		//vector<thread> threads;
		for (int i = 0; i < starts.size(); ++i)
		{
			int start = starts[i];
			int partition = vertexmap[start];
			graphs[partition].bfs_search(start, i, start, 0);
			//threads.push_back(thread(&search_act, start, i, &graphs[partition]));
		}
		//for (int i = 0; i < threads.size(); ++i)
		//	threads[i].join();
	}
	int getmultival(int index, int i)
	{
		int partition = vertexmap[index];
		return graphs[partition].getval(index,i);
	}
	vector<SubGraph> graphs;
	vector<SubGraph*> neighbors;
	vector<int> vertexmap;
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

	int criticallength = 0;

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
				if (count > criticallength)
					criticallength = count;
				//cout << endl;
			}
		}
	}
	cout << "Critical length of " << criticallength << endl;
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
					sample = hw.getval(sample,j);
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
				  }
				  if (counts[j][i].second != count)
					  cout << "Counts not equal! " << count << " for bfs != " << counts[j][i].second << " for pbfs!\n";
				  //cout << endl;
			  }
		  }
	  }
		/*
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
					}
					if (counts[j][i].second != count)
						cout << "Counts not equal! " << count << " for pbfs != " << counts[j][i].second << " for bfs!\n";
					//cout << endl;
				}
			}
		}*/
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