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

typedef boost::property<first_name_t, std::tuple<int, bool, int, int>> Color; //parent, color, partition, distance
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
void bfs_search_async(Subgraph* s, int index, int parent, int dist);

typedef hpx::lcos::local::spinlock mutex_type;

struct pennant
{
	pennant*left;
	pennant* right;
	int value;
	int dist;
	int parent;
	pennant()
	{
		value = 0;
		left = nullptr;
		right = nullptr;
	};
	pennant(int val)
	{
		value = val;
		left = nullptr;
		right = nullptr;
	};

};

pennant* merge(pennant* x, pennant* y)
{
	y->right = x->left;
	x->left = y;
	return y;
}
pennant* split(pennant* x)
{
	pennant* y = x->left;
	x->left = y->right;
	y->right = nullptr;
	return y;
}



typedef pair<pennant*, pennant*> fapair;

#define falsepen nullptr

fapair FA(pennant* x, pennant* y, pennant* z)
{
	if (x != nullptr && y == nullptr && z == nullptr)
	{
		return fapair(x, falsepen);
	}
	if (x == nullptr && y != nullptr && z == nullptr)
	{
		return fapair(y, falsepen);
	}
	if (x == nullptr && y == nullptr && z != nullptr)
	{
		return fapair(z, falsepen);
	}
	if (x != nullptr && y != nullptr && z == nullptr)
	{
		return fapair(falsepen, merge(x, y));
	}
	if (x != nullptr && y == nullptr && z != nullptr)
	{
		return fapair(falsepen, merge(x, z));
	}
	if (x == nullptr && y != nullptr && z != nullptr)
	{
		return fapair(falsepen, merge(y, z));
	}
	if (x != nullptr && y != nullptr && z != nullptr)
	{
		return fapair(x, merge(y, z));
	}
	//false,false,false
	return fapair(falsepen, falsepen);
}
struct bag
{
	vector<pennant*> values;
	int range;
	bag(int allocate)
	{
		values = vector<pennant*>(allocate, falsepen);
		range = allocate;
	};
	bag()
	{
		range = 0;
	};
	bool empty()
	{
		for (int i = 0; i < values.size(); ++i)
		{
			if (values[i] != nullptr)
				return false;
		}
		return true;
	}
	int size()
	{
		int count = 0;
		for (int i = 0; i < values.size(); ++i)
		{
			if (values[i] != nullptr)
				++count;
		}
		return count;
	}
	void insert(pennant* x)
	{
		int k = 0;
		while (values[k] != nullptr)
		{
			x = merge(values[k], x);
			values[k++] = nullptr;
		}
		values[k] = x;
	}
	void bag_merge(bag b)
	{
		pennant* y = falsepen;
		for (int k = 0; k < range; ++k)
		{
			fapair p = FA(values[k], b.values[k],y);
			values[k] = p.first;
			y = p.second;
		}
	}
	bag bag_split()
	{
		bag s = bag(range);
		pennant* y = values[0];
		values[0] = falsepen;
		for (int k = 1; k < range; ++k)
		{
			if (values[k] != nullptr)
			{
				s.values[k - 1] = values[k];
				values[k - 1] = values[k];
				values[k] = falsepen;
			}
		}
		if (y != nullptr)
			insert(y);
		return s;
	}
};


struct Subgraph
{
	Graph g;
	std::vector<pennant> pennants;
	int grainsize = 9999;
	mutex_type* m;
	void bfs_search(int index)
	{
		property_map < Graph, vertex_index_t >::type
			index_map = get(vertex_index, g);
		pennants = std::vector <pennant>(num_vertices(g));
		for (int i = 0; i < pennants.size(); ++i)
		{
			pennants[i].value = i;
			pennants[i].dist = 999999999;
			pennants[i].parent = -1;
		}
		pennants[index].parent = index;
		pennants[index].dist = 0;
		std::queue<int> q;
		q.push(index);
		int dist = 0;

		while (!q.empty())
		{
			int ind = q.back();
			q.pop();
			++dist;
			int parent = ind;
			if (pennants[ind].dist < dist)
				dist = pennants[ind].dist;
			graph_traits < Graph >::adjacency_iterator ai, a_end;

			for (boost::tie(ai, a_end) = adjacent_vertices(ind, g); ai != a_end; ++ai)
			{
				int ind = get(index_map, *ai);
				{
					if (pennants[ind].dist > dist + 1)
					{
						pennants[ind].dist = dist + 1;
						pennants[ind].parent = parent;
						q.push(ind);
					}
				}
			}
					
		}
	};

	void pbfs_search(int index)
	{
		pennants = std::vector <pennant> (num_vertices(g));
		for (int i = 0; i < pennants.size(); ++i)
		{
			pennants[i].value = i;
			pennants[i].dist = 999999999;
			pennants[i].parent = -1;
		}
		pennants[index].parent = index;
		pennants[index].dist = 0;
		vector<bag> v;
		int dist = 0;
		v.push_back(bag(33));
		v.back().insert(&pennants[index]);
		while (!v[dist].empty())
		{
			v.push_back(bag(33));
			process_layor(v[dist], v[dist + 1], dist);
			++dist;
		}
	};
	void process_layor(bag& in_bag, bag& out_bag, int d)
	{
		{
			property_map < Graph, vertex_index_t >::type
				index_map = get(vertex_index, g);
			for (int i = in_bag.values.size()-1; i >= 0; --i)
			{
				if (in_bag.values[i] != nullptr)
				{
					int parent = in_bag.values[i]->value;
					int dist = d;
					graph_traits < Graph >::adjacency_iterator ai, a_end;
					boost::tie(ai, a_end) = adjacent_vertices(parent, g);

					for (; ai != a_end; ++ai)
					{
						int ind = get(index_map, *ai);
						if (pennants[ind].dist > 99999)
						{
							//mutex_type::scoped_lock l(*m);
							pennants[ind].dist = dist + 1;
							pennants[ind].parent = parent;
							out_bag.insert(&(pennants[ind]));
						}
					}
				}
			}
			return;
		}
		return;
	};

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
  int start = randnodes(rng);
  Subgraph sub;
  sub.g = g;
  sub.grainsize = 32;
  sub.pbfs_search(start);
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
		  sample = sub.pennants[sample].parent;
	  }
	  //cout << endl;
  }
  sub.bfs_search(start);
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
		  sample = sub.pennants[sample].parent;
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