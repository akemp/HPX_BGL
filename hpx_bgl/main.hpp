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

using namespace boost;
using namespace std;
struct GraphComponent :
    hpx::components::managed_component_base<GraphComponent>
{

    static void parreset(MultiGraph* g, int start, int size, int toggled)
    {
        property_map<MultiGraph, multi_name_t>::type
            name = get(multi_name_t(), *g);
        for (int i = start; i < start + size; ++i)
        {
            name[i] = vector<int>(toggled, -1);
        }
    }

    void reset(int toggled)
    {
        int adder = grainsize;

        vector<hpx::future<void>> futs;
        MultiGraph* ptr = &g;

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
    void set(vector<vector<int>> nodes, int size, int edge, int starts)
    {
        for (int i = 0; i < nodes.size(); ++i)
        {
            for (int j = 0; j < nodes[i].size(); ++j)
                add_edge(i, nodes[i][j], g);
        }
        grainsize = size;
        edgefact = edge;

        name = get(multi_name_t(), g);
        multireset(starts);
    }
    int getval(int i, int index)
    {
        return name[i][index];
    }
    void multireset(int starts)
    {
        for (int i = 0; i < num_vertices(g); ++i)
        {
            name[i] = vector<int>(starts, -1);
        }
    }

    void setmulti(vector<vector<int>> nodes, int size, int edge, int starts)
    {
        for (int i = 0; i < nodes.size(); ++i)
        {
            for (int j = 0; j < nodes[i].size(); ++j)
                add_edge(i, nodes[i][j], g);
        }
        grainsize = size;
        edgefact = edge;

        name = get(multi_name_t(), g);
        multireset(starts);
    }

    void multival(vector<int> starts, bool sequential)
    {
        if (!sequential)
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
        else
        {
            for (int i = 0; i < starts.size(); ++i)
            {
                mpbfs(starts[i], i);
            }
        }
    }
    static void runMp(vector<int>::iterator loc, int index, int size, GraphComponent* gc)
    {
        for (int i = 0; i < size; ++i)
        {
            gc->mpbfs(*loc, index + i);
            ++loc;
        }
    }
    static vector<int> process_layor_multi(int loc, vector<int> in_bag, MultiGraph* g)
    {
        property_map < MultiGraph, vertex_index_t >::type
            index_map = get(vertex_index, *g);
        property_map<MultiGraph, multi_name_t>::type
            name = get(multi_name_t(), *g);
        vector<int> out_bag;
        int count = 0;
        for (int i = 0; i < in_bag.size(); ++i)
        {
            int val = in_bag[i];
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
    void mpbfs(int index, int loc)
    {
        name[index][loc] = index;
        vector<int> v;
        int dist = 0;
        v.push_back(index);
        MultiGraph* ptr = &g;
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
                        futures.push_back(hpx::async(hpx::util::bind(&process_layor_multi, loc, vector<int>(it, it + grainsize), ptr)));
                    else
                    {
                        futures.push_back(hpx::async(hpx::util::bind(&process_layor_multi, loc, vector<int>(it, it + (v.size() - last)), ptr)));
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
    int getmultival(int index, int i)
    {
        return name[index][i];
    }

    void bfs_search_act(vector<int> starts)
    {

        for (int i = 0; i < starts.size(); ++i)
        {
            bfs_search(starts[i], i);
        }
    }

    void bfs_search(int index, int loc)
    {
        property_map < MultiGraph, vertex_index_t >::type
            index_map = get(vertex_index, g);
        //pennants = std::vector <int>(num_vertices(g), -1);
        name[index][loc] = index;
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
            graph_traits < MultiGraph >::adjacency_iterator ai, a_end;

            for (boost::tie(ai, a_end) = adjacent_vertices(ind, g); ai != a_end; ++ai)
            {
                int ind = get(index_map, *ai);
                if (name[ind][loc] < 0)
                {
                    name[ind][loc] = parent;
                    q.push_back(ind);
                }
            }

        }
    };

    int getnum()
    {
        return num_vertices(g);
    }

    HPX_DEFINE_COMPONENT_ACTION(GraphComponent, set, set_action);
    HPX_DEFINE_COMPONENT_ACTION(GraphComponent, getval, getval_action);

    HPX_DEFINE_COMPONENT_ACTION(GraphComponent, setmulti, setmulti_action);
    HPX_DEFINE_COMPONENT_ACTION(GraphComponent, multival, multival_action);
    HPX_DEFINE_COMPONENT_ACTION(GraphComponent, getmultival, getmultival_action);

    HPX_DEFINE_COMPONENT_ACTION(GraphComponent, bfs_search_act, bfs_search_action);

    HPX_DEFINE_COMPONENT_ACTION(GraphComponent, getnum, num_vertices_action);
    property_map<MultiGraph, multi_name_t>::type
        name;
    MultiGraph g;
    int grainsize = 9999;
    int edgefact = 16;
    bool active = false;
};

typedef managed_component<GraphComponent> server_type;
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(server_type, GraphComponent);

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

typedef GraphComponent::bfs_search_action bfs_search_action;
HPX_REGISTER_ACTION_DECLARATION(bfs_search_action);
HPX_REGISTER_ACTION(bfs_search_action);

typedef GraphComponent::num_vertices_action num_vertices_action;
HPX_REGISTER_ACTION_DECLARATION(num_vertices_action);
HPX_REGISTER_ACTION(num_vertices_action);

struct graph_manager : client_base<graph_manager, GraphComponent>
{
    typedef client_base<graph_manager, GraphComponent> base_type;

    graph_manager(hpx::future<hpx::id_type> && id) : base_type(std::move(id)) {}

    void pbfs_search(vector<int >index)
    {
        hpx::async<multival_action>(this->get_gid(), index, true).get();
    }
    void set(vector<vector<int>>& edges, int grainsize, int edgefact, int starts)
    {
        hpx::async<set_action>(this->get_gid(), edges, grainsize, edgefact, starts).get();
    }
    int getval(int index, int i)
    {
        return hpx::async<getmultival_action>(this->get_gid(), index, i).get();
    }
    void setmulti(vector<vector<int>>& edges, int grainsize, int edgefact, int starts)
    {
        hpx::async<set_action>(this->get_gid(), edges, grainsize, edgefact, starts).get();
    }
    void multival(vector<int> starts)
    {
        hpx::async<multival_action>(this->get_gid(), starts, false).get();
    }

    void bfs_search(vector<int> starts)
    {
        hpx::async<bfs_search_action>(this->get_gid(), starts).get();
    }
    int getmultival(int index, int i)
    {

        return hpx::async<getmultival_action>(this->get_gid(), index, i).get();
    }
    int getnum()
    {

        return hpx::async<num_vertices_action>(this->get_gid()).get();
    }

};