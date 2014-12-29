//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "headers.hpp"

#include <hpx/include/components.hpp>

using hpx::components::stub_base;
using hpx::components::client_base;
using hpx::components::managed_component;
using hpx::components::managed_component_base;

typedef std::pair < unsigned int, unsigned int > Edge;
typedef std::vector<Edge> Edges;
using namespace boost;
using namespace std;

struct first_name_t {
	typedef boost::vertex_property_tag kind;
};

typedef boost::property<first_name_t, int> Color; //parent, color, partition, distance
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	Color> Graph;
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

void createEdges(const std::vector<std::vector<idx_t>>& nodes, Edges& edges)
{

	for (int i = 0; i < nodes.size(); ++i)
	{
		for (int j = 0; j < nodes[i].size(); ++j)
		{
			edges.push_back(Edge(i, nodes[i][j]));
		}
	}
	return;
}

typedef struct Subgraph;

typedef hpx::lcos::local::spinlock mutex_type;

vector<int> process_layor(vector<int>::iterator in_bag, Graph* g, int grainsize);

struct GraphComponent : 
	hpx::components::managed_component_base<GraphComponent>
{
	
	void pbfs_search(int index)
	{
		//pennants = std::vector <int> (num_vertices(g), -1);
		name[index] = index;
		vector<int> v;
		int dist = 0;
		v.push_back(index);
		Graph* ptr = &g;
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
						futures.push_back(hpx::async(hpx::util::bind(&process_layor, it, ptr, grainsize)));
					else
					{
						futures.push_back(hpx::async(hpx::util::bind(&process_layor, it, ptr, v.size() - last)));
						break;
					}
				}
			}
			//hpx::wait_all(futures.begin(), futures.end());
			vector<int> children;
			children.reserve(v.size() * edgefact);
			for (int i = 0; i < futures.size(); ++i)
			{
				vector<int> future = futures[i].get();
				children.insert(children.end(), future.begin(), future.end());
			}
			v = children;
		}
	}

	void set(Edges edges, int size, int edge)
	{
		for (int i = 0; i < edges.size(); ++i)
			add_edge(edges[i].first, edges[i].second, g);
		grainsize = size;
		edgefact = edge;

		name = get(first_name_t(), g);
		for (int i = 0; i < num_vertices(g); ++i)
		{
			name[i] = -1;
		}
	}
	void reset()
	{
		for (int i = 0; i < num_vertices(g); ++i)
		{
			name[i] = -1;
		}
	}
	int getval(int index)
	{
		return name[index];
	}

	HPX_DEFINE_COMPONENT_ACTION(GraphComponent, pbfs_search, pbfs_search_action);
	HPX_DEFINE_COMPONENT_ACTION(GraphComponent, set, set_action);
	HPX_DEFINE_COMPONENT_ACTION(GraphComponent, getval, getval_action);
	HPX_DEFINE_COMPONENT_ACTION(GraphComponent, reset, reset_action);
public:
	property_map<Graph, first_name_t>::type
		name;
private:
	Graph g;
	int grainsize = 9999;
	int edgefact = 16;
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


typedef GraphComponent::reset_action reset_action;
HPX_REGISTER_ACTION_DECLARATION(reset_action);
HPX_REGISTER_ACTION(reset_action);


struct graph_manager : client_base<graph_manager, GraphComponent>
{
	typedef client_base<graph_manager, GraphComponent> base_type;

	graph_manager(hpx::future<hpx::id_type> && id) : base_type(std::move(id)) {}

	void pbfs_search(int index)
	{
		hpx::async<pbfs_search_action>(this->get_gid(),index).get();
	}
	void set(Edges edges, int grainsize, int edgefact)
	{
		hpx::async<set_action>(this->get_gid(), edges, grainsize, edgefact).get();
	}
	int getval(int index)
	{
		return hpx::async<getval_action>(this->get_gid(), index).get();
	}
	void reset()
	{
		hpx::async<reset_action>(this->get_gid()).get();
	}
	
};


struct Subgraph
{
	Graph g;
	property_map<Graph, first_name_t>::type
		name;
	//std::vector<int> pennants;
	int grainsize = 9999;
	int edgefact = 16;
	void reset()
	{
		name = get(first_name_t(), g);
		for (int i = 0; i < num_vertices(g); ++i)
			name[i] = -1;
	}
	void bfs_search(int index)
	{
		property_map < Graph, vertex_index_t >::type
			index_map = get(vertex_index, g);
		//pennants = std::vector <int>(num_vertices(g), -1);
		name[index] = index;
		vector<int> q;
		q.reserve(num_vertices(g));
		q.push_back(index);
		int dist = 0;
		int spot = 0;
		while (spot < q.size())
		{
			int ind = q[spot];
			++spot;
			int parent = ind;
			graph_traits < Graph >::adjacency_iterator ai, a_end;

			for (boost::tie(ai, a_end) = adjacent_vertices(ind, g); ai != a_end; ++ai)
			{
				int ind = get(index_map, *ai);
				if (name[ind] < 0)
				{
					name[ind] = parent;
					q.push_back(ind);
				}
			}
					
		}
	};

	void pbfs_search(int index)
	{
		//pennants = std::vector <int> (num_vertices(g), -1);
		name[index] = index;
		vector<int> v;
		int dist = 0;
		v.push_back(index);
		Graph* ptr = &g;
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
						futures.push_back(hpx::async(hpx::util::bind(&process_layor, it, ptr, grainsize)));
					else
					{
						futures.push_back(hpx::async(hpx::util::bind(&process_layor, it, ptr, v.size() - last)));
						break;
					}
				}
			}
			//hpx::wait_all(futures.begin(), futures.end());
			vector<int> children;
			children.reserve(v.size() * edgefact);
			for (int i = 0; i < futures.size(); ++i)
			{
				vector<int> future = futures[i].get();
				children.insert(children.end(), future.begin(), future.end());
			}
			v = children;
		}
	};

};

vector<int> process_layor(vector<int>::iterator in_bag, Graph* g, int grainsize)
{
	property_map < Graph, vertex_index_t >::type
		index_map = get(vertex_index, *g);
	property_map<Graph, first_name_t>::type
		name = get(first_name_t(), *g);
	vector<int> out_bag;
	int count = 0;
	for (int i = 0; i < grainsize; ++i)
	{
		int val = *in_bag;
		++in_bag;
		graph_traits < Graph >::adjacency_iterator ai, a_end;

		for (boost::tie(ai, a_end) = adjacent_vertices(val, *g); ai != a_end; ++ai)
		{
			int ind = get(index_map, *ai);
			if (name[ind] >= 0)
				continue;
			name[ind] = val;
			out_bag.push_back(ind);
		}
	}
	return out_bag;
};
///////////////////////////////////////////////////////////////////////////////
int main()
{
	
	using namespace std;
	cout << "threads: " << hpx::get_os_thread_count() << endl;

	int nnodes;
	cout << "Enter nodes: ";
	cin >> nnodes;
	int ind;
	cout << "Enter edges per node: ";
	cin >> ind;

	vector<vector<idx_t>> nodes(nnodes);
	boost::random::mt19937 rng;
	boost::random::uniform_int_distribution<> randnodes(0, nodes.size()-1);
	for (int i = 0; i < nodes.size(); ++i)
	{
		vector<idx_t> edger;
		for (int j = 0; j < ind; ++j)
		{
			int spot = randnodes(rng);
			if (spot != i && std::find(nodes[spot].begin(), nodes[spot].end(), i) == nodes[spot].end())
				edger.push_back(spot);
			else
				--j;
		}
		nodes[i] = edger;
	}
	Edges edges;
	createEdges(nodes, edges);

	int grainsize = 128;

  Graph g(nodes.size());
  for (int i = 0; i < edges.size(); ++i)
	  add_edge(edges[i].first, edges[i].second, g);
  int start = randnodes(rng);
  Subgraph sub;
  sub.g = g;
  sub.grainsize = grainsize;
  sub.edgefact = ind;
  sub.reset();
  hpx::util::high_resolution_timer t;
  t.restart();
  sub.pbfs_search(start);
  double elapsed = t.elapsed();
  cout << elapsed << "s for parallel\n";
  vector<pair<int,int>> counts(32,pair<int,int>(-1, 0));
  for (int i = 0; i < 32; ++i)
  {
	  int sample = randnodes(rng);
	  counts[i].first = sample;
	  while (sample != start)
	  {
		  counts[i].second++;
		  if (sample == -1)
		  {
			  counts[i].second = -1;
			  break;
		  }
		  //cout << sample << "-" << sub.pennants[sample].dist << " ";
		  sample = sub.name[sample];
	  }
	  //cout << endl;
  }

  graph_manager hw = graph_manager::create(hpx::find_here());
  hw.set(edges, grainsize, ind);
  t.restart();
  hw.pbfs_search(start);
  elapsed = t.elapsed();
  cout << elapsed << "s for parallel component\n";

  for (int i = 0; i < 32; ++i)
  {
	  int sample = counts[i].first;
	  int count = 0;
	  while (sample != start)
	  {
		  count++;
		  if (sample == -1)
		  {
			  count = -1;
			  break;
		  }
		  //cout << sample << "-" << sub.pennants[sample].dist << " ";
		  sample = hw.getval(sample);
	  }
	  if (counts[i].second != count)
		  cout << "Counts not equal! " << count << " for bfs != " << counts[i].second << " for pbfs!\n";
	  //cout << endl;
  }

  sub.reset();
  t.restart();
  sub.bfs_search(start);
  elapsed = t.elapsed();
  cout << elapsed << "s for serial\n";
  for (int i = 0; i < 32; ++i)
  {
	  int sample = counts[i].first;
	  int count = 0;
	  while (sample != start)
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
	  if (counts[i].second != count)
		  cout << "Counts not equal! " << count << " for bfs != " << counts[i].second << " for pbfs!\n";
	  //cout << endl;
  }


  int s;
  cin >> s;
	return 0;
}