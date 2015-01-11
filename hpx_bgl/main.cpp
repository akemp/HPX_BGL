#include "main.hpp"

void parallel_edge_gen(vector<packed_edge>::iterator pedges, vector<vector<int>>* nodes, int size, vector<mutex_type*>* muts)
{

    boost::random::mt19937 rng;
    boost::random::uniform_int_distribution<> randnodes(0, nodes->size()-1);
    rng.seed(pedges->v0_low);
	for (int i = 0; i < size; ++i)
	{
        int v0 = randnodes(rng);//pedges->v0_low;
        int v1 = randnodes(rng);//pedges->v1_low;
        ++pedges;
		if (v0 == v1)
			continue;
		{
			//undirected so no changes to final edgelist
			if (v1 < v0)
			{
                int temp = v0;
				v0 = v1;
				v1 = temp;
			}
			mutex_type::scoped_lock l1(*(*muts)[v0]);
			if (std::find((*nodes)[v0].begin(), (*nodes)[v0].end(), v1) != (*nodes)[v0].end())
				continue;
			(*nodes)[v0].push_back(v1);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
int main()
{

	using namespace std;
	cout << "threads: " << hpx::get_os_thread_count() << endl;

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
	boost::random::mt19937 rng;
	boost::random::uniform_int_distribution<> randnodes(0, nnodes*edgefactor);

	uint64_t seed1 = 2, seed2 = 3;
    uint_fast32_t seed[5];

	make_mrg_seed(seed1, seed2, seed);

	//Edges edges;
	vector<vector<int>> nodes(
		nnodes
		);
	vector<mutex_type*> muts(nodes.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		muts[i] = new mutex_type;
	}
    uint64_t ind = (uint64_t)(edgefactor) << scale;
	vector<packed_edge> pedges(ind);
    generate_kronecker_range(seed, scale, 0, (uint64_t)pedges.size(), &pedges.front());
	cout << "Kronecker range generated. Making edgelist.\n";
	{

		vector<hpx::thread> edgefuts;
		int addsize = grainsize * 16;
		for (int i = 0; i < pedges.size(); i += addsize)
		{
			if (i + addsize < pedges.size())
			{
				edgefuts.push_back(hpx::thread(
					hpx::util::bind(&parallel_edge_gen, pedges.begin() + i, &nodes, addsize, &muts))
					);
			}
			else
			{
				edgefuts.push_back(hpx::thread(
					hpx::util::bind(&parallel_edge_gen, pedges.begin() + i, &nodes, pedges.size() - i, &muts))
					);
				break;
			}
		}
		for (int i = 0; i < edgefuts.size(); ++i)
			edgefuts[i].join();
	}
	cout << "Edgelist generated. Running tests.\n";


	vector<int> starts(searches);


	vector<vector<pair<int, int>>> counts;
	if (acctest != 0)
	{
		cout << "Setting up serial values for testing.\n";
		counts = vector<vector<pair<int, int>>>(starts.size(), vector<pair<int, int>>(64, pair<int, int>(-1, 0)));
	}
	int nverts;
	{
		graph_manager hw = graph_manager::create(hpx::find_here());
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
			//sub.reset();
			//sub.bfs_search(starts[j]);
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
		

		{
			hpx::util::high_resolution_timer t;
			t.restart();
		  graph_manager hw = graph_manager::create(hpx::find_here());
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
	  }
	  {
		  hpx::util::high_resolution_timer t;
		  t.restart();
		  graph_manager hw = graph_manager::create(hpx::find_here());
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

		cout << "Accuracy tests complete.\n";
	}
	cout << "Serial comparison:\n";
	{
		hpx::util::high_resolution_timer t;
		t.restart();
		graph_manager hw = graph_manager::create(hpx::find_here());
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