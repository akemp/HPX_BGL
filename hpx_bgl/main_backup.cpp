//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "headers.hpp"

#include <hpx/include/components.hpp>

#include "../generator/make_graph.h"
#include "../generator/utils.h"

using hpx::components::stub_base;
using hpx::components::client_base;
using hpx::components::managed_component;
using hpx::components::managed_component_base;

typedef hpx::lcos::local::spinlock mutex_type;

struct multi_name_t {
	typedef boost::vertex_property_tag kind;
};

typedef boost::property<multi_name_t, vector<int> > MultiColor; //parent, color, partition, distance
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	MultiColor> MultiGraph;

typedef multi_name_t lock_name_t;
typedef MultiGraph LockGraph;

int get_parts(std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, std::vector<idx_t> &part,
	idx_t nparts)
{
	// Needed by parmetis
	idx_t nvtxs = xadj.size() - 1, ncon = 1;
	idx_t *vsize = NULL;
	idx_t *vwgt = NULL, *adjwgt = NULL;
	real_t *tpwgts = NULL, *ubvec = NULL;
	idx_t *options = NULL;
	idx_t objval;
	int result = METIS_PartGraphKway(
		&nvtxs, &ncon, &xadj[0], &adjncy[0], vwgt,
		vsize, adjwgt, &nparts, tpwgts,
		ubvec, options, &objval, &part[0]);
	return result;
}
void toCSR(const std::vector<std::vector<idx_t>>& nodes, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy)
{
	int total = 0;
	xadj.push_back(0);
	for (int i = 0; i < nodes.size(); ++i)
	{
		for (int j = 0; j < nodes[i].size(); ++j)
		{
			adjncy.push_back(nodes[i][j]);
			++total;
		}
		xadj.push_back(total);
	}
	return;
}

typedef struct Subgraph;

struct GraphComponent :
	hpx::components::managed_component_base<GraphComponent>
{
	static void parreset(LockGraph* g, int start, int size, int toggled)
	{
		property_map<LockGraph, lock_name_t>::type
			name = get(lock_name_t(), *g);
		for (int i = start; i < start + size; ++i)
		{
			name[i] = vector<int>(toggled, -1);
		}
	}
	void pbfs_search(int index, int toggled)
	{
		name[index][toggled] = index;
		vector<int> v;
		int dist = 0;
		v.push_back(index);
		LockGraph* ptr = &g;
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
						futures.push_back(hpx::async(&process_layor_component, it, ptr, grainsize, toggled));
					else
					{
						futures.push_back(hpx::async(&process_layor_component, it, ptr, v.size() - last, toggled));
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

	void reset(int toggled)
	{
		int adder = grainsize;

		vector<hpx::future<void>> futs;
		LockGraph* ptr = &g;

		for (int i = 0; i < num_vertices(g); i += adder)
		{
			if (i + adder < num_vertices(g))
			{
				futs.push_back(hpx::async(&parreset, ptr, i, adder, toggled));
			}
			else
			{
				futs.push_back(hpx::async(&parreset, ptr, i, num_vertices(g) - i, toggled));
			}
		}
		hpx::wait_all(futs);

	}
	void set(const vector<vector<uint32_t>> nodes, int size, int edge, int nverts, int toggled)
	{
		g = LockGraph(nverts);
		grainsize = size;
		edgefact = edge;
		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int j = 0; j < nodes[i].size(); ++j)
				add_edge(i, nodes[i][j], g);
		}

		name = get(lock_name_t(), g);
		reset(toggled);
	}
	int getval(int i, int index)
	{
		return name[index][i];
	}
	void multireset(int starts)
	{
		for (int i = 0; i < num_vertices(gs); ++i)
		{
			names[i] = vector<int>(starts, -1);
		}
	}

	void setmulti(vector<vector<uint32_t>> nodes, int size, int edge, int starts)
	{
		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int j = 0; j < nodes[i].size(); ++j)
				add_edge(i, nodes[i][j], gs);
		}
		grainsize = size;
		edgefact = edge;

		names = get(multi_name_t(), gs);
		multireset(starts);
	}

	void multival(vector<int> starts)
	{
		vector<hpx::future<void>> futures;
		int i = 0;
		int adder = 1;
		for (vector<int>::iterator it = starts.begin(); it < starts.end(); it += adder)
		{
			int last = i;
			i += adder;
			if (i < starts.size())
			{
				futures.push_back(hpx::async(&runMp, it, last, adder, this));
			}
			else
			{
				futures.push_back(hpx::async(&runMp, it, last, starts.size() - last, this));
				break;
			}
		}
		hpx::wait_all(futures);
	}
	static void runMp(vector<int>::iterator loc, int index, int size, GraphComponent* gc)
	{
		for (int i = 0; i < size; ++i)
		{
			gc->mpbfs(index + i, *loc);
			++loc;
		}
	}
	static vector<int> process_layor_component(vector<int>::iterator in_bag,
		LockGraph* g, const int grainsize, const int toggle)
	{
		vector<int> out_bag;


		const property_map < LockGraph, vertex_index_t >::type
			index_map = get(vertex_index, *g);
		const property_map<LockGraph, lock_name_t>::type
			name = get(lock_name_t(), *g);
		for (int i = 0; i < grainsize; ++i)
		{
			int val = *in_bag;
			++in_bag;
			graph_traits < LockGraph >::adjacency_iterator ai, a_end;
			boost::tie(ai, a_end) = adjacent_vertices(val, *g);
			for (; ai != a_end; ++ai)
			{
				int ind = get(index_map, *ai);
				if (name[toggle][ind] >= 0)
					continue;
				name[toggle][ind] = val;
				out_bag.push_back(ind);
			}
		}
		return out_bag;
	}

	static vector<int> process_layor_multi(int loc, vector<int>::iterator in_bag, MultiGraph* g, int grainsize)
	{
		property_map < MultiGraph, vertex_index_t >::type
			index_map = get(vertex_index, *g);
		property_map<MultiGraph, multi_name_t>::type
			name = get(multi_name_t(), *g);
		vector<int> out_bag;
		int count = 0;
		for (int i = 0; i < grainsize; ++i)
		{
			int val = *in_bag;
			++in_bag;
			graph_traits < MultiGraph >::adjacency_iterator ai, a_end;

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
	void mpbfs(int loc, int index)
	{
		//pennants = std::vector <int> (num_vertices(g), -1);
		names[index][loc] = index;
		vector<int> v;
		int dist = 0;
		v.push_back(index);
		MultiGraph* ptr = &gs;
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
						futures.push_back(hpx::async(hpx::util::bind(&process_layor_multi, loc, it, ptr, grainsize)));
					else
					{
						futures.push_back(hpx::async(hpx::util::bind(&process_layor_multi, loc, it, ptr, v.size() - last)));
						break;
					}
				}
			}
			//hpx::wait_all(futures.begin(), futures.end());
			vector<int> children;
			for (int i = 0; i < futures.size(); ++i)
			{
				vector<int> future = futures[i].get();
				children.insert(children.end(), future.begin(), future.end());
			}
			v = children;
		}
	}
	int getmultival(int i, int index)
	{
		return names[index][i];
	}

	HPX_DEFINE_COMPONENT_ACTION(GraphComponent, pbfs_search, pbfs_search_action);
	HPX_DEFINE_COMPONENT_ACTION(GraphComponent, set, set_action);
	HPX_DEFINE_COMPONENT_ACTION(GraphComponent, getval, getval_action);

	HPX_DEFINE_COMPONENT_ACTION(GraphComponent, setmulti, setmulti_action);
	HPX_DEFINE_COMPONENT_ACTION(GraphComponent, multival, multival_action);
	HPX_DEFINE_COMPONENT_ACTION(GraphComponent, getmultival, getmultival_action);
public:
	property_map<LockGraph, lock_name_t>::type
		name;
	property_map<MultiGraph, multi_name_t>::type
		names;
private:
	LockGraph g;
	MultiGraph gs;
	int grainsize = 9999;
	int edgefact = 16;
	bool active = false;
};


typedef managed_component<GraphComponent> server_type;
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(server_type, GraphComponent);

typedef GraphComponent::pbfs_search_action pbfs_search_action;
HPX_REGISTER_ACTION_DECLARATION(pbfs_search_action);
HPX_REGISTER_ACTION(pbfs_search_action);

typedef GraphComponent::set_action set_action;
HPX_REGISTER_ACTION_DECLARATION(set_action);
HPX_REGISTER_ACTION(set_action);

typedef GraphComponent::getval_action getval_action;
HPX_REGISTER_ACTION_DECLARATION(getval_action);
HPX_REGISTER_ACTION(getval_action);

typedef GraphComponent::setmulti_action setmulti_action;
HPX_REGISTER_ACTION_DECLARATION(setmulti_action);
HPX_REGISTER_ACTION(setmulti_action);


typedef GraphComponent::multival_action multival_action;
HPX_REGISTER_ACTION_DECLARATION(multival_action);
HPX_REGISTER_ACTION(multival_action);

typedef GraphComponent::getmultival_action getmultival_action;
HPX_REGISTER_ACTION_DECLARATION(getmultival_action);
HPX_REGISTER_ACTION(getmultival_action);


struct graph_manager : client_base<graph_manager, GraphComponent>
{
	typedef client_base<graph_manager, GraphComponent> base_type;

	graph_manager(hpx::future<hpx::id_type> && id) : base_type(std::move(id)) {}

	void pbfs_search(int index, int toggle)
	{
		hpx::async<pbfs_search_action>(this->get_gid(), index, toggle).get();
	}
	void set(const vector<vector<uint32_t>>& nodes, int grainsize, int edgefact, int nverts, int toggled)
	{
		//	void set(const vector<vector<uint32_t>> nodes, int size, int edge, int nverts, int toggled)
		hpx::async<set_action>(this->get_gid(), nodes, grainsize, edgefact, nverts, toggled).get();
	}
	int getval(int index, int i)
	{
		return hpx::async<getval_action>(this->get_gid(), i, index).get();
	}
	void setmulti(vector<vector<uint32_t>>& edges, int grainsize, int edgefact, int starts)
	{
		hpx::async<setmulti_action>(this->get_gid(), edges, grainsize, edgefact, starts).get();
	}
	void multival(vector<int> starts)
	{

		hpx::async<multival_action>(this->get_gid(), starts).get();
	}
	int getmultival(int index, int i)
	{

		return hpx::async<getmultival_action>(this->get_gid(), i, index).get();
	}

};


void parallel_edge_gen(vector<packed_edge>::iterator pedges, vector<vector<uint32_t>>* nodes, int size, vector<mutex_type*>* muts)
{
	for (int i = 0; i < size; ++i)
	{
		int v0 = pedges->v0_low;
		int v1 = pedges->v1_low;
		++pedges;
		if (v0 == v1)
			continue;
		{
			//undirected so no changes to final edgelist
			if (v1 > v0)
			{
				int temp = v0;
				v0 = v1;
				v1 = v0;
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

	int nnodes;
	int scale;
	cout << "Enter scale (nodes=2^input): ";
	cin >> scale;
	int ind;
	cout << "Enter edges per node: ";
	cin >> ind;
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
	int acctest;
	cout << "Run accuracy tests (0 for no, 1 for yes)?";
	cin >> acctest;
	nnodes = pow(2, scale);
	//
	boost::random::mt19937 rng;
	boost::random::uniform_int_distribution<> randnodes(0, nnodes*ind);

	uint64_t seed1 = 4, seed2 = 3;

	uint_fast32_t seed;

	make_mrg_seed(seed1, seed2, &seed);

	//Edges edges;
	vector<vector<uint32_t>> nodes(
		nnodes
		);
	vector<mutex_type*> muts(nodes.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		muts[i] = new mutex_type;
	}
	vector<packed_edge> pedges(nnodes*ind);
	generate_kronecker_range(&seed, scale, 0, pedges.size(), &pedges.front());
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
		Subgraph sub;
		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int j = 0; j < nodes[i].size(); ++j)
				add_edge(i, nodes[i][j], sub.g);
		}
		nverts = num_vertices(sub.g);
		randnodes = boost::random::uniform_int_distribution<>(0, nverts - 1);
		sub.grainsize = grainsize;
		sub.edgefact = ind;
		for (int j = 0; j < starts.size(); ++j)
		{
			starts[j] = randnodes(rng);
			if (acctest != 0)
			{
				sub.reset();
				sub.bfs_search(starts[j]);
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
							break;
						}
						//cout << sample << "-" << sub.pennants[sample].dist << " ";
						sample = sub.name[sample];
					}
					counts[j][i].second = count;
					//cout << endl;
				}
			}
		}
	}
	if (acctest != 0)
	{
		cout << "Running accuracy tests.\n";
		{
			Subgraph sub;
			for (int i = 0; i < nodes.size(); ++i)
			{
				for (int j = 0; j < nodes[i].size(); ++j)
					add_edge(i, nodes[i][j], sub.g);
			}

			sub.grainsize = grainsize;
			sub.edgefact = ind;
			for (int j = 0; j < starts.size(); ++j)
			{
				sub.reset();
				sub.pbfs_search(starts[j]);
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
							sample = sub.name[sample];
						}
						if (counts[j][i].second != count)
							cout << "Counts not equal! " << count << " for bfs != " << counts[j][i].second << " for pbfs!\n";
						//cout << endl;
					}
				}
			}
		}

  {
	  graph_manager hw = graph_manager::create(hpx::find_here());
	  hw.set(nodes, grainsize, ind, nverts, starts.size());
	  for (int j = 0; j < starts.size(); ++j)
	  {
		  hw.pbfs_search(starts[j], j);
		  
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
				  sample = hw.getval(j,sample);
			  }
			  if (counts[j][i].second != count)
				  cout << "Counts not equal! " << count << " for bfs != " << counts[j][i].second << " for pbfs!\n";
			  //cout << endl;
		  }

	  }
  }
  {
	  graph_manager hw = graph_manager::create(hpx::find_here());
	  hw.setmulti(nodes, grainsize, ind, starts.size());
	  hw.multival(starts);

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
				  sample = hw.getmultival(j, sample);
			  }
			  if (counts[j][i].second != count)
				  cout << "Counts not equal! " << count << " for bfs != " << counts[j][i].second << " for pbfs!\n";
			  //cout << endl;
		  }
	  }
  }

		cout << "Accuracy tests complete.\n";
	}
	cout << "Benchmarking.\n";
	{
		hpx::util::high_resolution_timer t;
		t.restart();
		Subgraph sub;
		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int j = 0; j < nodes[i].size(); ++j)
				add_edge(i, nodes[i][j], sub.g);
		}


		sub.grainsize = grainsize;
		sub.edgefact = ind;
		hpx::util::high_resolution_timer t1;
		for (int j = 0; j < starts.size(); ++j)
		{
			sub.reset();
			sub.bfs_search(starts[j]);
		}
		double elapsed = t.elapsed();
		double elapsed1 = t1.elapsed();
		cout << elapsed << "s for serial\n";
		cout << elapsed1 << "s search time for serial\n";
	}

  {
	  hpx::util::high_resolution_timer t;
	  t.restart();
	  Subgraph sub;
	  for (int i = 0; i < nodes.size(); ++i)
	  {
		  for (int j = 0; j < nodes[i].size(); ++j)
			  add_edge(i, nodes[i][j], sub.g);
	  }

	  sub.grainsize = grainsize;
	  sub.edgefact = ind;
	  hpx::util::high_resolution_timer t1;
	  for (int j = 0; j < starts.size(); ++j)
	  {
		  sub.reset();
		  sub.pbfs_search(starts[j]);
	  }
	  double elapsed = t.elapsed();
	  double elapsed1 = t1.elapsed();
	  cout << elapsed << "s for parallel\n";
	  cout << elapsed1 << "s search time for parallel\n";
  }
  {

	  hpx::util::high_resolution_timer t;
	  t.restart();
	  graph_manager hw = graph_manager::create(hpx::find_here());
	  hw.set(nodes, grainsize, ind, nverts, starts.size());
	  hpx::util::high_resolution_timer t1;
	  for (int j = 0; j < starts.size(); ++j)
	  {
		  hw.pbfs_search(starts[j], j);
	  }
	  double elapsed = t.elapsed();
	  double elapsed1 = t1.elapsed();
	  cout << elapsed << "s for parallel component\n";
	  cout << elapsed1 << "s search time for parallel component\n";
  }
  {
	  hpx::util::high_resolution_timer t;
	  t.restart();
	  graph_manager hw = graph_manager::create(hpx::find_here());
	  hw.setmulti(nodes, grainsize, ind, starts.size());
	  hpx::util::high_resolution_timer t1;
	  hw.multival(starts);
	  double elapsed = t.elapsed();
	  double elapsed1 = t1.elapsed();
	  cout << elapsed << "s for highly parallel component\n";
	  cout << elapsed1 << "s search time for highly parallel component\n";
  }

	int s;
	cin >> s;
	return 0;
}