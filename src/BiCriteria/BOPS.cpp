#include <memory>
#include <algorithm>
#include <time.h>

#include "BOPS.h"

#define SCREEN false

BOPS::BOPS(const AdjacencyMatrix &adj_matrix, Pair<double> eps, Heuristic &h_f, Heuristic &h_b,
    const LoggerPtr logger, size_t perimeter_factor) : 
    BOAStar(adj_matrix, eps, logger), heuristic_f(h_f), heuristic_b(h_b) {
    set_perimeter_factor(perimeter_factor);
}

// TODO: only count the 1st element of the front to front heuristic, and count the 2nd if needed
void BOPS::operator() (size_t source, size_t target, Heuristic &heuristic, SolutionSet & solutions,
    unsigned int time_limit) {

    // Set the perimeter as the h1 from the target to the source / factor
    vector<size_t> tar_h = heuristic_b(target);
    if (SCREEN) cout << "heuristic_b(target)[0]" << tar_h.front() << endl;
    set_perimeter_dist(tar_h[0]/perimeter_factor);

    // Gcl
    bool global_dominance = false, local_dominance = false;
    vector<vector_heap> g_cl;
    for (size_t i = 0; i < adj_matrix.size()+1; i++) {
        vector_heap tmp_g_close;
        g_cl.push_back(tmp_g_close);
    }

    start_time = clock();
    this->start_logging(source, target);

    NodePtr node_f;
    NodePtr node_b;
    NodePtr next;

    // Saving all the unused NodePtrs in a vector improves performace for some reason
    vector<NodePtr> closed_f;  // forward
    vector<NodePtr> closed_b;  // backward

    // Vector to hold mininum cost of 2nd criteria per node
    vector<size_t> min_g2_f(this->adj_matrix.size()+1, MAX_COST);
    vector<size_t> min_g2_b(this->adj_matrix.size()+1, MAX_COST);

    // Init open heap
    Node::more_than_full_cost more_than;
    vector<NodePtr> open_f;
    make_heap(open_f.begin(), open_f.end(), more_than);
    vector<NodePtr> open_b;
    make_heap(open_b.begin(), open_b.end(), more_than);

    node_f = make_shared<Node>(source, vector<size_t>(2,0), heuristic_f(source), nullptr, true);
    node_f->set_path();

    node_b = make_shared<Node>(target, vector<size_t>(2,0), vector<size_t>(2,0), nullptr, false);
    node_b->set_path();

    node_f->set_h_node(node_b);
    node_b->set_h_node(node_f);

    open_f.push_back(node_f);
    std::push_heap(open_f.begin(), open_f.end(), more_than);
    list<PathGvalPair> open_f_paths = get_paths(open_f);

    open_b.push_back(node_b);
    std::push_heap(open_b.begin(), open_b.end(), more_than);
    list<PathGvalPair> open_b_paths = get_paths(open_b);

    // Generate the perimeter
    while (!open_b.empty()) {
        if ((clock() - start_time)/CLOCKS_PER_SEC > time_limit) {
            this->end_logging(solutions, false);
            return;
        }

        // Pop min from queue and process
        std::pop_heap(open_b.begin(), open_b.end(), more_than);
        node_b = open_b.back();
        open_b.pop_back();

        if (SCREEN) {
            assert(node_b->other_h.empty());
            cout << "---------------------------------" << endl;
            cout << "expand node_b" << endl;
            cout << *node_b << endl;
            cout << "min_g2_b[source]: " << min_g2_b[source] << endl;
            cout << "min_g2_b[node_b->id]: " << min_g2_b[node_b->id] << endl;
            cout << "---------------------------------" << endl;
        }

        // Dominance check (only local)
        // If we are going to use sets of heuristics, then we need to change checking rule
        if (node_b->g[1] >= min_g2_b[node_b->id]) {
            reinsert(node_b, open_b, closed_b, more_than);
            continue;
        }

        // We decide to expand this node
        // TODO: might need the other heuristics for the pathmax heuristic 
        // if front to front is too time-consuming
        min_g2_b[node_b->id] = node_b->g[1];
        num_expansion_b ++;
        num_expansion ++;

        closed_b.push_back(node_b);  // We fully expand this node, so add to closed list

        if (node_b->id == source) {
            // Find one solution during the search (not put into the perimeter)
            solutions.push_back(node_b);
            num_sol_b ++;
        } else if (node_b->g[0] < perimeter_dist) {
            // If g1 of the node_b is smaller than the perimeter_dist, then expand
            // Check to which neighbors we should extend the paths
            const vector<Edge> &outgoing_edges = adj_matrix[node_b->id];
            branching_factor_b += outgoing_edges.size() - 1;
            if (SCREEN)
                cout << "outgoing_edges size: " << outgoing_edges.size() << endl;
            for (auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                // Ignore the node we are coming from
                if (node_b->parent != nullptr && p_edge->target == node_b->parent->id) continue;

                size_t next_id = p_edge->target;
                vector<size_t> next_g = node_b->g + p_edge->cost;
                if (next_g[1] >= min_g2_b[next_id]) continue;  // Dominance check

                next = make_shared<Node>(next_id, next_g, vector<size_t>(next_g.size(), 0), 
                    node_b, node_b->is_forward);
                next->set_path();
                open_b.push_back(next);
                push_heap(open_b.begin(), open_b.end(), more_than);
                num_generation ++;
                num_generation_b ++;

                if (SCREEN) {
                    cout << "\tgenerate node" << endl;
                    cout << "\t" << *next << endl;
                }
            }
        } else {
            // Add node_b to the perimeter
            // TODO: Might conside how to sort the node in perimeter
            perimeter.push_back(node_b);
        }
        if (SCREEN) {
            cout << "iteration done, open_b.size: " << open_b.size() << ", closed_b.size: "
                << closed_b.size() << endl;
            cout << endl;
        }
    }
    runtime_pre_h += clock() - start_time;  // End the backward search

    // Begin with forward search
    if (SCREEN)
    {
        cout << "Start forward search with perimeter size: " << perimeter.size() << endl;
        cout << "Num of sol from backward: " << solutions.size() << endl;
        cout << "Update the heuristics for the source" << endl;
    }

    // Update the heuristic for the backward open list given the forward open list
    clock_t tmp_start = clock();
    for (auto& _node_ : open_f) {
        list<HeuristicNodePair> tmp_h_list;
        for (const auto& _p_node_ : perimeter) {
            assert(!_p_node_->is_cand);
            vector<size_t> tmp_h = get_diff_heuristic(_node_->id, _p_node_->id);  // front 2 front
            tmp_h += _p_node_->g;
            bool ignore_h = false;  // Check if the heuristic is dominanted by the heuristic
            list<HeuristicNodePair>::iterator h_it = tmp_h_list.begin();
            while (h_it != tmp_h_list.end()) {
                if (is_bounded(h_it->first, tmp_h)) {
                    h_it = tmp_h_list.erase(h_it);
                } else if (is_bounded(tmp_h, h_it->first)) {
                    ignore_h = true;
                    break;
                } else {
                    ++ h_it;
                }
            }

            if (!ignore_h) {
                tmp_h_list.emplace_back(tmp_h, _p_node_);
            }
        }
        _node_->set_h_val(tmp_h_list);
    }
    runtime_update_h += clock() - tmp_start;

    while (!open_f.empty()) {
        if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit) {
            this->end_logging(solutions, false);
            return;
        }

        // Pop min from queue and process
        pop_heap(open_f.begin(), open_f.end(), more_than);
        node_f = open_f.back();
        open_f.pop_back();

        if (SCREEN) {
            cout << "Current node: " << *node_f << endl;
        }

        // Global dominance check, TODO: binary search
        assert(!global_dominance);
        for (const auto& sol : solutions) {
            assert(sol->g[0] == sol->f[0] && sol->g[1] == sol->f[1]);
            if (is_bounded(node_f->f, sol->f)) {
                if (SCREEN) 
                    cout << "\tglobal dominant: " << *node_f << endl;
                reinsert(node_f, open_f, closed_f, more_than);
                global_dominance = true;
                break;
            } else if (node_f->f[0] < sol->f[0]) {
                break;
            }
        }
        if (global_dominance) {
            global_dominance = false;
            continue;
        }

        // Local dominance check, TODO: binary search
        assert(!local_dominance);
        if (SCREEN) print_gcl(g_cl, node_f->id);
        for (const auto& g_min_val : g_cl[node_f->id]) {
            if (is_bounded(node_f->g, g_min_val)) {
                if (SCREEN) cout << "\tlocal  dominant: " << *node_f << endl;
                reinsert(node_f, open_f, closed_f, more_than);
                local_dominance = true;
                break;
            } else if (node_f->g[0] < g_min_val[0]) {
                break;
            }
        }
        if (local_dominance) {
            local_dominance = false;
            continue;
        }

        // If the current node reaches the perimeter, then
        // 1. Add to the solution set
        // 2. Switch to the next heuristic in other_h
        // 3. Do not update g_cl since we do not officially expand this node
        // 4. continue to the next iteration
        if (node_f->id == node_f->h_node->id) {
            // This is a solution, will be checked in the opposite search
            // Don't need to put in the node to the OPEN
            NodePtr sol = make_shared<Node>(node_f->id, node_f->g + node_f->h_node->g, 
                vector<size_t>(node_f->g.size(), 0), node_f->parent, true, node_f->h_node, true);
            sol->combine_path();
            solutions.push_back(sol);
            num_sol_f ++;
            // min_f2 = sol->g[1];

            if (SCREEN) {
                cout << "\tThis is a solution" << endl;
                cout << "\t" << *sol << endl;
                cout << "\t----- open_f -----" << endl;
                print_list(open_f, more_than, 5, true);
                cout << "\t--- end open_f ---" << endl;
            }

            // Reinsert the node with another heuristic to the open
            reinsert(node_f, open_f, closed_f, more_than);
            continue;
        }

        // We decide to expand this node
        // TODO: might need the other heuristics for the pathmax heuristic 
        // if front to front is too time-consuming
        if (SCREEN) {
            cout << "\tExpand" << endl;
            cout << "\t----- open_f -----" << endl;
            print_list(open_f, more_than, 5, true);
            cout << "\t--- end open_f ---" << endl;
            cout << endl;
        }
        g_cl[node_f->id].push(node_f->g);
        // if (g_cl[node_f->id].size() > 1)
        // {
        //     cout << "\t";
        //     print_gcl(g_cl, node_f->id);
        //     cout << endl;
        // } 
        num_expansion ++;
        num_expansion_f ++;

        // Check to which neighbors we should extend the paths
        const vector<Edge> &outgoing_edges = adj_matrix[node_f->id];
        branching_factor_f += outgoing_edges.size() - 1;
        if (SCREEN)
            cout << "\t\toutgoing_edges size: " << outgoing_edges.size() << endl;
        for (auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            // Ignore the node we are coming from
            if (node_f->parent != nullptr && p_edge->target == node_f->parent->id) continue;
            size_t next_id = p_edge->target;
            vector<size_t> next_g = node_f->g + p_edge->cost;

            // Local dominance check
            assert(!local_dominance);
            for (const auto& g_min_val : g_cl[next_id]) {
                if (is_bounded(next_g, g_min_val)) {
                    local_dominance = true;
                    break;
                } else if (next_g[0] < g_min_val[0]) {
                    break;
                }
            }
            if (local_dominance) {
                local_dominance = false;
                continue;
            }

            // Update the heuristic for the child node
            clock_t update_h_start = clock();
            list<HeuristicNodePair> h_list;
            for (const auto& other_node : perimeter) {
                vector<size_t> next_h = get_diff_heuristic(next_id, other_node->id);  // f2f

                // // Debug: check if f2f heuristic is admissible
                // if (SCREEN && next_id != source) {
                //     // Generate inverse graph
                //     vector<Edge> tmp_edges;
                //     for (size_t v = 0; v < adj_matrix.get_graph_size(); v++) {
                //         vector<Edge> tmp_tmp_edges = adj_matrix[v+1];
                //         for (const auto& tmp_tmp_e : tmp_tmp_edges) {
                //             tmp_edges.push_back(tmp_tmp_e);
                //         }
                //     }
                //     AdjacencyMatrix inv_graph(adj_matrix.get_graph_size(), tmp_edges, true);
                //     ShortestPathHeuristic sp_heuristic(next_id, adj_matrix.get_graph_size(), inv_graph);
                //     using std::placeholders::_1;
                //     Heuristic tmp_tmp_h = bind(&ShortestPathHeuristic::operator(), sp_heuristic, _1);
                //     vector<size_t> tmp_tmp_h1 = tmp_tmp_h(other_node->id);

                //     cout << "next_h: " << next_h[0] << ", " << next_h[1] << endl;
                //     cout << "tmp_tmp_h1: " << tmp_tmp_h1[0] << ", " << tmp_tmp_h1[1] << endl;
                //     cout << endl;
                //     assert(next_h[0] <= tmp_tmp_h1[0]);
                //     assert(next_h[1] <= tmp_tmp_h1[1]);
                // }
                // // End debug

                next_h += other_node->g;
                vector<size_t> next_f = next_g + next_h;

                // // Dominance check
                // if ((next_g[1]+next_h[1]) >= min_f2 || (next_g[1] >= min_g2_f[next_id])) {
                //     continue;
                // }

                // Global dominance check
                assert(!global_dominance);
                for (const auto& sol : solutions) {
                    assert(sol->g[0] == sol->f[0] && sol->g[1] == sol->f[1]);
                    if (is_bounded(next_f, sol->f)) {
                        global_dominance = true;
                        break;
                    } else if (next_f[0] < sol->f[0]) {
                        break;
                    }
                }
                if (global_dominance) {
                    global_dominance = false;
                    continue;
                }

                // Check if the heuristic is dominant, either the current heuristic dominant 
                // one of the existing hueristics or the other way around
                bool ignore_h = false;
                list<HeuristicNodePair>::iterator h_it = h_list.begin();
                while (h_it != h_list.end()) {
                    if (is_bounded(h_it->first, next_h)) {
                        h_it = h_list.erase(h_it);
                    } else if (is_bounded(next_h, h_it->first)) {
                        ignore_h = true;
                        break;
                    } else {
                        ++ h_it;
                    }
                }
                if (ignore_h) continue;
                h_list.emplace_back(next_h, other_node);
            }
            runtime_update_h += clock() - update_h_start;

            if (!h_list.empty()) {
                // If not dominated create node and push to queue
                // Creation is defered after dominance check as it is
                // relatively computational heavy and should be avoided if possible
                clock_t node_gen_start = clock();
                next = make_shared<Node>(next_id, next_g, h_list, node_f, node_f->is_forward);
                next->set_path();
                open_f.push_back(next);
                push_heap(open_f.begin(), open_f.end(), more_than);

                // Statistic analysis
                num_generation ++;
                num_generation_f ++;
                runtime_node_gen += clock() - node_gen_start;
            
                if (SCREEN) {
                    cout << "\tgenerate node" << endl;
                    cout << "\t" << *next << endl;
                    // next->print_h();
                }
            }
        }
        closed_f.push_back(node_f);  // We fully expand this node, so add to closed list
        if (SCREEN) {
            cout << "\t\tgeneration done, open_f.size: " << open_f.size() << ", closed_f.size: "
                << closed_f.size() << endl;
            cout << endl;
        }
    }
}

vector<size_t> BOPS::get_diff_heuristic(size_t loc1, size_t loc2) {
    if (loc1 == loc2) {
        size_t num_of_obj = heuristic_f(loc1).size();
        return vector<size_t>(num_of_obj, 0);
    }
    vector<size_t> h1 = heuristic_f(loc1);
    vector<size_t> h2 = heuristic_f(loc2);
    vector<size_t> tmp_h1 = get_diff(h1, h2);

    h1 = heuristic_b(loc1);
    h2 = heuristic_b(loc2);
    vector<size_t> tmp_h2 = get_diff(h1, h2);
    return get_comax(tmp_h1, tmp_h2);
}

void BOPS::reinsert(NodePtr node, vector<NodePtr>& open, vector<NodePtr>& closed,
    const Node::more_than_full_cost& more_than) {
    if (!node->other_h.empty()) {
        node->update_h();
        open.push_back(node);
        push_heap(open.begin(), open.end(), more_than);
        if (SCREEN)
        {
            cout << "\tReinsert: " << *node << " to OPEN" << endl;
            cout << "\t--- open list after reinseert ---" << endl;
            print_list(open, more_than, 5, true);
            cout << "\t--------- end reinseert ---------" << endl;
        }
    } else {
        closed.push_back(node);
        if (SCREEN)
        {
            cout << "Reinsert: " << *node << " to CLOSED" << endl;
        }
    }
}

void BOPS::print_gcl(const vector<vector_heap>& in_gcl, int id) {
    cout << "g_cl[" << id << "]: ";
    
    vector_heap tmp_gcl(in_gcl[id]);
    while (!tmp_gcl.empty()) {
        vector<size_t> top_g = tmp_gcl.top();
        tmp_gcl.pop();
        cout << "(";
        for (size_t i = 0; i < top_g.size(); i++) {
            cout << top_g[i];
            if (i < top_g.size()-1) {
                cout << ", ";
            }
        }
        cout << ") ";
    }
    cout << endl;
}