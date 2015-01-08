//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include "headers.hpp"
#include "main.hpp"

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
	{
		cout << "Running accuracy tests.\n";
		{

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