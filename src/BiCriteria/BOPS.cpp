#include <memory>
#include <algorithm>
#include <time.h>

#include "BOPS.h"

BOPS::BOPS(const AdjacencyMatrix &adj_matrix, Pair<double> eps, Heuristic &h_f, Heuristic &h_b,
    double perimeter_factor, uint b_mode, const LoggerPtr logger, uint screen) : 
    BOAStar(adj_matrix, eps, logger, screen), heuristic_f(h_f), 
    heuristic_b(h_b), backward_mode(b_mode) {
    set_perimeter_factor(perimeter_factor);
}

// TODO: only count the 1st element of the front to front heuristic, and count the 2nd if needed
void BOPS::operator() (size_t source, size_t target, Heuristic &heuristic, SolutionSet & solutions,
    uint time_limit) {

    // Set the perimeter as the h1 from the target to the source / factor
    vector<size_t> tar_h = heuristic_b(target);
    set_perimeter_dist((double) tar_h[0] * perimeter_factor);
    if (screen > DEBUG_LOG_EXPANSION) 
        cout << "perimeter_dist: " << perimeter_dist << endl;

    // Gcl
    bool global_dominance = false, local_dominance = false;
    vector<vector_heap> g_cl;
    for (size_t i = 0; i < adj_matrix.size()+1; i++) {
        vector_heap tmp_g_close;
        g_cl.push_back(tmp_g_close);
    }

    clock_t tmp_start = clock();
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

    if (backward_mode == 0)
        node_b = make_shared<Node>(target, vector<size_t>(2,0), vector<size_t>(2,0), nullptr, false);
    else
        node_b = make_shared<Node>(target, vector<size_t>(2,0), heuristic_b(target), nullptr, false);
    node_b->set_path();

    node_f->set_h_node(node_b);
    node_b->set_h_node(node_f);

    open_f.push_back(node_f);
    push_heap(open_f.begin(), open_f.end(), more_than);
    list<PathGvalPair> open_f_paths = get_paths(open_f);

    open_b.push_back(node_b);
    push_heap(open_b.begin(), open_b.end(), more_than);
    list<PathGvalPair> open_b_paths = get_paths(open_b);

    // Generate the perimeter
    if (screen > DEBUG_LOG_EXPANSION) cout << "Generate perimeter... ";
    while (!open_b.empty()) {
        if ((clock() - start_time)/CLOCKS_PER_SEC > time_limit) {
            if (screen > DEBUG_LOG_EXPANSION) this->end_logging(solutions, false);
            return;
        }

        // Pop min from queue and process
        pop_heap(open_b.begin(), open_b.end(), more_than);
        node_b = open_b.back();
        open_b.pop_back();

        if (screen > DEBUG_LOG_DETAILES) {
            assert(node_b->other_h.empty());
            cout << "\n---------------------------------" << endl;
            cout << "expand node_b ";
            cout << *node_b << endl;
            if (min_g2_b[source] < SIZE_MAX) 
                cout << "\tmin_g2_b[source]: " << min_g2_b[source] << endl;
            if (min_g2_b[node_b->id] < SIZE_MAX) 
                cout << "\tmin_g2_b[" << node_b->id << "]: " << min_g2_b[node_b->id] << endl;
            cout << "---------------------------------" << endl;
        }

        // Dominance check
        if ((((1+this->eps[1])*node_b->f[1]) >= min_g2_b[source]) ||
            (node_b->g[1] >= min_g2_b[node_b->id])) {
            closed_b.push_back(node_b);
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
        } else if ((double) node_b->f[0] < perimeter_dist) {
            // If f1 (or g1) of the node_b is smaller than the perimeter_dist, then expand
            // Check to which neighbors we should extend the paths
            const vector<Edge> &outgoing_edges = adj_matrix[node_b->id];
            branching_factor_b += outgoing_edges.size() - 1;
            if (screen > DEBUG_LOG_DETAILES)
                cout << "outgoing_edges size: " << outgoing_edges.size() << endl;
            for (auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                // Ignore the node we are coming from
                if (node_b->parent != nullptr && p_edge->target == node_b->parent->id) continue;

                size_t next_id = p_edge->target;
                vector<size_t> next_g = node_b->g + p_edge->cost;
                vector<size_t> next_h;
                switch (backward_mode) {
                case 0:  // BOD
                    next_h = vector<size_t>(next_g.size(), 0);
                    break;
                case 1: // BOA
                    next_h = heuristic_b(source);
                    break;
                default:
                    exit(-1);
                    break;
                }

                if ((next_g[1]+next_h[1]) >= min_g2_b[source] || next_g[1] >= min_g2_b[next_id]) {
                    continue;
                }

                next = make_shared<Node>(next_id, next_g, next_h, node_b, false);
                next->set_h_node(node_f);
                next->set_path();
                open_b.push_back(next);
                push_heap(open_b.begin(), open_b.end(), more_than);
                num_generation ++;
                num_generation_b ++;

                if (screen > DEBUG_LOG_DETAILES) {
                    cout << "\tgenerate node" << endl;
                    cout << "\t" << *next << endl;
                }
            }
        } else {
            // Add node_b to the perimeter
            // TODO: Might conside how to sort the node in perimeter
            perimeter.push_back(node_b);
        }
    }
    runtime_pre_h += clock() - start_time;  // End the backward search

    // Begin with forward search
    if (screen > DEBUG_LOG_EXPANSION)
    {
        cout << "done, open_b.size: " << open_b.size() << ", closed_b.size: "
            << closed_b.size() << ", #sol: " << solutions.size() 
            << ", perimeter size: " << perimeter.size() << ", f1: " << perimeter[0]->f[0] << endl;
        cout << "Start forward search... " << endl;
        cout << "Update the heuristics for the source... " << endl;
    }

    // Update the heuristic for the backward open list given the forward open list
    tmp_start = clock();
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
    if (screen > DEBUG_LOG_EXPANSION) cout << "done" << endl;

    while (!open_f.empty()) {
        if ((clock() - start_time)/CLOCKS_PER_SEC > time_limit) {
            if (screen > DEBUG_LOG_EXPANSION) this->end_logging(solutions, false);
            return;
        }

        // Pop min from queue and process
        pop_heap(open_f.begin(), open_f.end(), more_than);
        node_f = open_f.back();
        open_f.pop_back();

        if (screen > DEBUG_LOG_EXPANSION) {
            cout << "Current node: " << *node_f;
        }

        // Global dominance check, TODO: binary search
        tmp_start = clock();
        assert(!global_dominance);
        for (const auto& sol : solutions) {
            assert(sol->g[0] == sol->f[0] && sol->g[1] == sol->f[1]);
            if (is_bounded(node_f->f, sol->f)) {
                reinsert(node_f, open_f, closed_f, more_than);
                global_dominance = true;
                break;
            } else if (node_f->f[0] < sol->f[0]) {
                break;
            }
        }
        runtime_global_check += clock() - tmp_start;
        if (global_dominance) {
            global_dominance = false;
            continue;
        }

        // Local dominance check, TODO: binary search
        tmp_start = clock();
        assert(!local_dominance);
        if (screen > DEBUG_LOG_DETAILES) print_gcl(g_cl, node_f->id);
        for (const auto& g_min_val : g_cl[node_f->id]) {
            if (is_bounded(node_f->g, g_min_val)) {
                reinsert(node_f, open_f, closed_f, more_than);
                local_dominance = true;
                break;
            } else if (node_f->g[0] < g_min_val[0]) {
                break;
            }
        }
        runtime_local_check += clock() - tmp_start;
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
            // sol->combine_path();
            solutions.push_back(sol);
            num_sol_f ++;

            if (screen > DEBUG_LOG_EXPANSION) {
                cout << "is a solution" << endl;
                cout << "\t" << *sol << endl;
            }

            // Reinsert the node with another heuristic to the open
            reinsert(node_f, open_f, closed_f, more_than);
            continue;
        }

        // We decide to expand this node
        // TODO: might need the other heuristics for the pathmax heuristic 
        // if front to front is too time-consuming
        if (screen > DEBUG_LOG_EXPANSION) cout << " is expanded" << endl;
        g_cl[node_f->id].push(node_f->g);
        num_expansion ++;
        num_expansion_f ++;

        // Check to which neighbors we should extend the paths
        const vector<Edge> &outgoing_edges = adj_matrix[node_f->id];
        branching_factor_f += outgoing_edges.size() - 1;
        if (screen > DEBUG_LOG_DETAILES)
            cout << "\toutgoing_edges size: " << outgoing_edges.size() << endl;
        for (auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            if ((clock() - start_time)/CLOCKS_PER_SEC > time_limit) {
                if (screen > DEBUG_LOG_EXPANSION) this->end_logging(solutions, false);
                return;
            }

            // Ignore the node we are coming from
            if (node_f->parent != nullptr && p_edge->target == node_f->parent->id) continue;
            clock_t node_gen_start = clock();
            size_t next_id = p_edge->target;
            vector<size_t> next_g = node_f->g + p_edge->cost;

            // Local dominance check
            tmp_start = clock();
            assert(!local_dominance);
            for (const auto& g_min_val : g_cl[next_id]) {
                if (is_bounded(next_g, g_min_val)) {
                    local_dominance = true;
                    break;
                } else if (next_g[0] < g_min_val[0]) {
                    break;
                }
            }
            runtime_local_check += clock() - tmp_start;
            if (local_dominance) {
                local_dominance = false;
                continue;
            }

            // Update the heuristic for the child node
            clock_t update_h_start = clock();
            list<HeuristicNodePair> h_list;
            for (const auto& other_node : perimeter) {
                vector<size_t> next_h = get_diff_heuristic(next_id, other_node->id);  // f2f
                next_h += other_node->g;
                vector<size_t> next_f = next_g + next_h;

                // Global dominance check
                tmp_start = clock();
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
                runtime_global_check += clock() - tmp_start;
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
                next = make_shared<Node>(next_id, next_g, h_list, node_f, node_f->is_forward);
                next->set_path();
                open_f.push_back(next);
                push_heap(open_f.begin(), open_f.end(), more_than);

                // Statistic analysis
                num_generation ++;
                num_generation_f ++;

                if (screen > DEBUG_LOG_EXPANSION) {
                    cout << "\tgenerate node " << *next << endl;
                    if (screen > DEBUG_LOG_DETAILES) next->print_h();
                }
            }
            runtime_node_gen += clock() - node_gen_start;
        }
        closed_f.push_back(node_f);  // We fully expand this node, so add to closed list
        if (screen > DEBUG_LOG_DETAILES) {
            cout << "\tgeneration done, open_f.size: " << open_f.size() << ", closed_f.size: "
                << closed_f.size() << endl;
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
        num_reinsert ++;
        if (screen > DEBUG_LOG_DETAILES)
        {
            cout << "\tReinsert: " << *node << " to OPEN" << endl;
            cout << "\t--- open list after reinseert ---" << endl;
            print_list(open, more_than, 5, true);
            cout << "\t--------- end reinseert ---------" << endl;
        }
    } else {
        closed.push_back(node);
        if (screen > DEBUG_LOG_DETAILES)
            cout << "Reinsert: " << *node << " to CLOSED" << endl;
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