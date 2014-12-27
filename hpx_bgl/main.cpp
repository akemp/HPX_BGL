//  Copyright (c) 2014 Andrew Kemp
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "headers.hpp"

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

struct Subgraph
{
	Graph g;
	//std::vector<int> pennants;
	int grainsize = 9999;
	void reset()
	{

		property_map<Graph, first_name_t>::type
			name = get(first_name_t(), g);
		for (int i = 0; i < num_vertices(g); ++i)
			name[i] = -1;
	}
	void bfs_search(int index)
	{
		property_map < Graph, vertex_index_t >::type
			index_map = get(vertex_index, g);
		property_map<Graph, first_name_t>::type
			name = get(first_name_t(), g);
		//pennants = std::vector <int>(num_vertices(g), -1);
		name[index] = index;
		list<int> q;
		q.push_back(index);
		int dist = 0;

		while (!q.empty())
		{
			int ind = q.front();
			q.pop_front();
			int parent = ind;
			graph_traits < Graph >::adjacency_iterator ai, a_end;

			for (boost::tie(ai, a_end) = adjacent_vertices(ind, g); ai != a_end; ++ai)
			{
				int ind = get(index_map, *ai);
				{
					if (name[ind] < 0)
					{
						name[ind] = parent;
						q.push_back(ind);
					}
				}
			}
					
		}
	};

	void pbfs_search(int index)
	{
		//pennants = std::vector <int> (num_vertices(g), -1);
		property_map<Graph, first_name_t>::type
			name = get(first_name_t(), g);
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
			children.reserve(v.size() * 16);
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

	vector<idx_t> xadj;
	vector<idx_t> adjncy;

	toCSR(nodes, xadj, adjncy);

	vector<idx_t> part(xadj.size() - 1);
	int partitions = 1;
	//cout << endl << "Enter partitions: ";
	//cin >> partitions;
	if (partitions > 1)
		get_parts(xadj, adjncy, part, partitions);
	else
	{
		partitions = 1;
		std::cout << "Serial execution.\n";
		std::fill(part.begin(), part.end(), 0);
	}

  Graph g(nodes.size());
  for (int i = 0; i < edges.size(); ++i)
	  add_edge(edges[i].first, edges[i].second, g);
  property_map<Graph, first_name_t>::type
	  name = get(first_name_t(), g);
  for (int i = 0; i < num_vertices(g); ++i)
	  name[i] = -1;
  int start = randnodes(rng);
  Subgraph sub;
  sub.g = g;
  sub.grainsize = 128;
  sub.reset();
  hpx::util::high_resolution_timer t;
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
		  sample = name[sample];
	  }
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
		  sample = name[sample];
	  }
	  if (counts[i].second != count)
		  cout << "Counts not equal! " << count << " for bfs != " << counts[i].second << " for pbfs!\n";
	  //cout << endl;
  }
  /*
  std::get<3>(name[start]) = 0;
  int loc = std::get<2>(name[start]);
  
  vector<Subgraph> subs(partitions);
  vector<Subgraph*> neighbors(partitions);
  vector<mutex_type> muts(partitions+1);
  for (int j = 0; j < partitions; ++j)
  {
	  neighbors[j] = &subs[j];//std::make_shared<Subgraph>(std::move(subs[j]));
  }
  for (int j = 0; j < partitions; ++j)
  {
	  subs[j].g = g;
	  subs[j].m = &muts[j];
  }
  hpx::util::high_resolution_timer t;
  subs[loc].bfs_search(start);
  double fin = t.elapsed();
  cout << "Completed in " << fin << "s\n";
  for (int j = 0; j < partitions; ++j)
  {

	  property_map<Graph, first_name_t>::type
		  tname = get(first_name_t(), subs[j].g);

	  for (int i = 0; i < num_vertices(g); ++i)
	  {
		  if (get<3>(name[i]) > get<3>(tname[i]) + 1)
		  {
			  get<0>(name[i]) = get<0>(tname[i]);
			  get<3>(name[i]) = get<3>(tname[i]) + 1;
		  }
	  }
  }
  for (int i = 0; i < 32; ++i)
  {
	  int sample = randnodes(rng);
	  while (sample != start)
	  {
		  if (sample == -1)
		  {
			  cout << "ERROR!";
			  break;
		  }
		  cout << sample << " ";
		  sample = get<0>(name[sample]);
	  }
	  cout << endl;
  }
  */
  int s;
  cin >> s;
	return 0;
}