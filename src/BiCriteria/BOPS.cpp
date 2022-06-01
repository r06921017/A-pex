#include <memory>
#include <algorithm>
#include <time.h>

#include "BOPS.h"

#define SCREEN false

BOPS::BOPS(const AdjacencyMatrix &adj_matrix, Pair<double> eps, Heuristic &h_f, Heuristic &h_b,
    const LoggerPtr logger, size_t lh_f, size_t lh_b) : 
    BOAStar(adj_matrix, eps, logger), lookahead_f(lh_f), lookahead_b(lh_b),
    heuristic_f(h_f), heuristic_b(h_b) {}

void BOPS::operator() (size_t source, size_t target, SolutionSet & solutions,
    unsigned int time_limit) {

    start_time = std::clock();
    this->start_logging(source, target);

    NodePtr node_f;
    NodePtr node_b;
    NodePtr next;
    list<NodePtr> tmp_path;

    // Store the paths
    list<PathGvalPair> closed_f_paths;
    list<PathGvalPair> closed_b_paths;

    // Saving all the unused NodePtrs in a vector improves performace for some reason
    std::vector<NodePtr> closed_f;  // forward
    std::vector<NodePtr> closed_b;  // backward

    // Vector to hold mininum cost of 2nd criteria per node
    std::vector<size_t> min_g2_f(this->adj_matrix.size()+1, MAX_COST);
    std::vector<size_t> min_g2_b(this->adj_matrix.size()+1, MAX_COST);

    // Init open heap
    Node::more_than_full_cost more_than;
    std::vector<NodePtr> open_f;
    std::make_heap(open_f.begin(), open_f.end(), more_than);
    std::vector<NodePtr> open_b;
    std::make_heap(open_b.begin(), open_b.end(), more_than);

    node_f = std::make_shared<Node>(source, std::vector<size_t>(2,0), heuristic_f(source), nullptr, true);
    node_f->set_path();
    open_f.push_back(node_f);
    std::push_heap(open_f.begin(), open_f.end(), more_than);
    list<PathGvalPair> open_f_paths = get_paths(open_f);

    node_b = std::make_shared<Node>(target, std::vector<size_t>(2,0), heuristic_b(target), nullptr, false);
    node_b->set_path();
    open_b.push_back(node_b);
    std::push_heap(open_b.begin(), open_b.end(), more_than);
    list<PathGvalPair> open_b_paths = get_paths(open_b);

    node_f->other_path = &open_b_paths.front();
    node_b->other_path = &open_f_paths.front();

    while (!open_f.empty() || !open_b.empty()) {
        // Begin with forward search
        cout << "Start forward search with open_f size: " << open_f.size() << endl;
        while (!open_f.empty()) {
            if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit) {
                this->end_logging(solutions, false);
                return;
            } else if (num_expansion_f > lookahead_f) {  // Switch to backward search
                num_expansion += num_expansion_f;
                num_generation += num_generation_f;
                num_expansion_f = 0;
                num_generation_f = 0;

                // extract paths from the current open list
                open_f_paths = get_paths(open_f);
                closed_f_paths = get_paths(closed_f);

                if (!closed_f_paths.empty()) {
                    // Remove path from forward OPEN if it is in the forward CLOSED
                    list<PathGvalPair>::iterator ofp_it = open_f_paths.begin();
                    while (ofp_it != open_f_paths.end()) {
                        if (find(closed_f_paths.begin(), closed_f_paths.end(), *ofp_it) != closed_f_paths.end()) {
                            ofp_it = open_f_paths.erase(ofp_it);
                        } else {
                            ++ ofp_it;
                        }
                    }
                }

                open_f.clear();
                closed_f.clear();
                update_open(open_b, solutions, target, open_b_paths, open_f_paths);

                if (SCREEN) {
                    cout << "\nForward search end" << endl;
                    cout << "open_b size: " << open_b.size() << endl;
                    cout << "open_b_paths size: " << open_b_paths.size() << endl;
                }
                break;
            }

            // Pop min from queue and process
            std::pop_heap(open_f.begin(), open_f.end(), more_than);
            node_f = open_f.back();
            open_f.pop_back();
            num_generation_f +=1;

            cout << "expand node_f " << num_generation_f << endl;

            if (SCREEN) {
                cout << "---------------------------------" << endl;
                cout << "expand node_f" << endl;
                cout << "id: " << node_f->id << endl;
                cout << "g: " << node_f->g[0] << ", " << node_f->g[1] << endl;
                cout << "h: " << node_f->h[0] << ", " << node_f->h[1] << endl;
                cout << "f: " << node_f->f[0] << ", " << node_f->f[1] << endl;
                cout << "path: ";
                for (const auto& p : node_f->path) {
                    cout << p << ", ";
                }
                cout << "\nother path: ";
                for (const auto& p : node_f->other_path->first) {
                    cout << p << ", ";
                }
                cout << "\nother g: " << node_f->other_path->second[0] << ", " 
                    << node_f->other_path->second[1] << endl;
                cout << "---------------------------------" << endl;
            }

            // Dominance check
            if ((((1+this->eps[1])*node_f->f[1]) >= min_g2_f[target]) ||
                (node_f->g[1] >= min_g2_f[node_f->id])) {
                closed_f.push_back(node_f);
                continue;
            }

            min_g2_f[node_f->id] = node_f->g[1];
            num_expansion_f += 1;

            // Find one solution during the search
            if (node_f->id == target) {
                // node_f->g = node_f->g + node_f->other_path->second;

                auto it = solutions.begin();
                while (it != solutions.end()) {
                    if (is_bounded(*it, node_f)) {
                        // candidate node is better than one of the solutions
                        it = solutions.erase(it);
                    } else {
                        ++it;
                    }
                }

                solutions.push_back(node_f);
                log_solution(node_f);

                cout << "find a solution from the forward: " << node_f->g[0] << ", " << node_f->g[1] << endl;
                num_sol_front ++;
                continue;
            }

            // Check to which neighbors we should extend the paths
            const std::vector<Edge> &outgoing_edges = adj_matrix[node_f->id];
            if (SCREEN)
                cout << "outgoing_edges size: " << outgoing_edges.size() << endl;
            for (auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                size_t next_id = p_edge->target;
                std::vector<size_t> next_g = node_f->g + p_edge->cost;

                vector<size_t> next_h = get_diff_heuristic(next_id, node_f->other_path->first.back());
                    next_h = next_h + node_f->other_path->second;

                // Dominance check
                if ((((1+this->eps[1])*(next_g[1]+next_h[1])) >= min_g2_f[target]) ||
                    (next_g[1] >= min_g2_f[next_id])) {
                    if (SCREEN)
                        cout << "is dominated in " << next_id << endl;    
                    continue;
                }

                // If not dominated create node and push to queue
                // Creation is defered after dominance check as it is
                // relatively computational heavy and should be avoided if possible
                next = std::make_shared<Node>(next_id, next_g, next_h, node_f,
                    node_f->is_forward, node_f->other_path);
                next->set_path();
                open_f.push_back(next);
                std::push_heap(open_f.begin(), open_f.end(), more_than);

                if (SCREEN) {
                    cout << "...................." << endl;
                    cout << "\tgenerate node" << endl;
                    cout << "\tid: " << next->id << endl;
                    cout << "\tg: " << next->g[0] << ", " << next->g[1] << endl;
                    cout << "\th: " << next->h[0] << ", " << next->h[1] << endl;
                    cout << "\tf: " << next->f[0] << ", " << next->f[1] << endl;
                    cout << "\tpath: ";
                    for (const auto& p : next->path) {
                        cout << p << ", ";
                    }
                    cout << "\n\tother path: ";
                    for (const auto& p : next->other_path->first) {
                        cout << p << ", ";
                    }
                    cout << "\n\tother g: " << next->other_path->second[0] << ", " 
                        << next->other_path->second[1] << endl;
                    cout << "...................." << endl;
                }
            }
            closed_f.push_back(node_f);
        }

        // Backward search
        // TODO: 1. extract the paths from OPEN, and then ignore the paths from CLOSED
        // 2. For each extracted path, update the front-to-front heuristic
        if (SCREEN)
            cout << "Start backward search with open_b size: " << open_b.size() << endl;
        while (!open_b.empty()) {
            if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit) {
                this->end_logging(solutions, false);
                return;
            } else if (num_expansion_b > lookahead_b) {  // Switch to forward search
                num_expansion += num_expansion_b;
                num_generation += num_generation_b;
                num_expansion_b = 0;
                num_generation_b = 0;

                // extract paths from the current open list
                open_b_paths = get_paths(open_b);
                closed_b_paths = get_paths(closed_b);

                if (!closed_b_paths.empty()) {
                    // Remove path from forward OPEN if it is in the forward CLOSED
                    list<PathGvalPair>::iterator obp_it = open_b_paths.begin();
                    while (obp_it != open_b_paths.end()) {
                        if (find(closed_b_paths.begin(), closed_b_paths.end(), *obp_it) != closed_b_paths.end()) {
                            obp_it = open_b_paths.erase(obp_it);
                        } else {
                            ++ obp_it;
                        }
                    }
                }
                open_b.clear();
                closed_b.clear();

                update_open(open_f, solutions, target, open_f_paths, open_b_paths);
                if (SCREEN) {
                    cout << "Backward search end" << endl;
                    cout << "open_f size: " << open_f.size() << endl;
                    cout << "open_f_paths size: " << open_f_paths.size() << endl;
                }
                break;
            }

            // Pop min from queue and process
            std::pop_heap(open_b.begin(), open_b.end(), more_than);
            node_b = open_b.back();
            open_b.pop_back();
            num_generation_b +=1;
            assert(node_b->path.back() == node_b->id);

            if (SCREEN) {
                cout << "---------------------------------" << endl;
                cout << "expand node_b" << endl;
                cout << "id: " << node_b->id << endl;
                cout << "g: " << node_b->g[0] << ", " << node_b->g[1] << endl;
                cout << "h: " << node_b->h[0] << ", " << node_b->h[1] << endl;
                cout << "f: " << node_b->f[0] << ", " << node_b->f[1] << endl;
                cout << "path: ";
                for (const auto& p : node_b->path) {
                    cout << p << ", ";
                }
                cout << "\nother path: ";
                for (const auto& p : node_b->other_path->first) {
                    cout << p << ", ";
                }
                cout << "\nother g: " << node_b->other_path->second[0] << ", " 
                    << node_b->other_path->second[1] << endl;
                cout << "---------------------------------" << endl;
            }

            // Dominance check
            if ((((1+this->eps[1])*node_b->f[1]) >= min_g2_b[source]) ||
                (node_b->g[1] >= min_g2_b[node_b->id])) {
                closed_b.push_back(node_b);
                continue;
            }

            min_g2_b[node_b->id] = node_b->g[1];
            num_expansion_b += 1;

            // Find one solution during the search
            if (node_b->id == source) {
                // node_b->g = node_b->g + node_b->other_path->second;

                auto it = solutions.begin();
                while (it != solutions.end()) {
                    if (is_bounded(*it, node_b)) {
                        // candidate node is better than one of the solutions
                        it = solutions.erase(it);
                    } else {
                        ++it;
                    }
                }

                solutions.push_back(node_b);
                log_solution(node_b);
                // if (SCREEN)
                    cout << "find a solution from the backward: " << node_b->g[0] << ", " << node_b->g[1] << endl;
                num_sol_back ++;
                continue;
            }

            // Check to which neighbors we should extend the paths
            const std::vector<Edge> &outgoing_edges = adj_matrix[node_b->id];
            if (SCREEN)
                cout << "outgoing_edges size: " << outgoing_edges.size() << endl;
            for (auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                size_t next_id = p_edge->target;
                std::vector<size_t> next_g = node_b->g + p_edge->cost;

                vector<size_t> next_h = get_diff_heuristic(next_id, node_b->other_path->first.back());
                    next_h = next_h + node_b->other_path->second;

                // Dominance check
                if ((((1+this->eps[1])*(next_g[1]+next_h[1])) >= min_g2_b[source]) ||
                    (next_g[1] >= min_g2_b[next_id])) {
                    if (SCREEN)
                        cout << "is dominated in " << next_id << endl;

                    continue;
                }

                // If not dominated create node and push to queue
                // Creation is defered after dominance check as it is
                // relatively computational heavy and should be avoided if possible
                next = std::make_shared<Node>(next_id, next_g, next_h, node_b, 
                    node_b->is_forward, node_b->other_path);
                next->set_path();
                open_b.push_back(next);
                std::push_heap(open_b.begin(), open_b.end(), more_than);

                if (SCREEN) {
                    cout << "...................." << endl;
                    cout << "\tgenerate node" << endl;
                    cout << "\tid: " << next->id << endl;
                    cout << "\tg: " << next->g[0] << ", " << next->g[1] << endl;
                    cout << "\th: " << next->h[0] << ", " << next->h[1] << endl;
                    cout << "\tf: " << next->f[0] << ", " << next->f[1] << endl;
                    cout << "\tpath: ";
                    for (const auto& p : next->path) {
                        cout << p << ", ";
                    }
                    cout << "\n\tother path: ";
                    for (const auto& p : next->other_path->first) {
                        cout << p << ", ";
                    }
                    cout << "\n\tother g: " << next->other_path->second[0] << ", " 
                        << next->other_path->second[1] << endl;
                    cout << "...................." << endl;
                }
            }
            closed_b.push_back(node_b);
        }
    }
}

// void BOPS::update_open(vector<NodePtr>& open, const vector<NodePtr>& other_open,
//     size_t target, SolutionSet &solutions) {

//     // Initialize the unodered_map that maps the path
//     vector<NodePtr> new_open;
//     map<list<size_t>, bool> paths_in_open;
//     for (const auto& tmp_node : open) {
//         if (paths_in_open.count(tmp_node->path) == 0) {
//             paths_in_open[tmp_node->path] = false;
//         }
//     }

//     for (const auto& node : open) {
//         // Check if we already update this path in the previous iteration
//         assert(paths_in_open.count(node->path) > 0);
//         if (paths_in_open[node->path]) {
//             continue;
//         } else {
//             paths_in_open[node->path] = true;
//         }

//         for (const auto& other_node : other_open) {
//             if (node->id == other_node->id) {  
//                 // This is a solution
//                 assert(node->path.back() == node->id);
//                 assert(other_node->path.back() == other_node->id);
//                 assert(node->id == other_node->id);

//                 // only for evaluation
//                 NodePtr cand_sol = std::make_shared<Node>(node->id, (node->g + other_node->g), 0, 
//                     node->parent, node->is_forward, other_node->parent);
//                 bool need_to_add = true;
//                 vector<NodePtr>::iterator it = solutions.begin();
//                 while (it != solutions.end()) {
//                     if (is_dominated_dr(cand_sol, *it)) {
//                         // candidate node is better than one of the solutions
//                         // TODO: check if this condition ever exists
//                         num_better_sol ++;
//                         it = solutions.erase(it);
//                     } else if (is_dominated_dr(*it, cand_sol)) {
//                         // one of the solutions is better than candidate node
//                         num_redundent_sol ++;
//                         need_to_add = false;
//                         break;
//                     } else {
//                         // The cadidate solution is not dominant by the current solution set
//                         ++it;
//                     }
//                 }
//                 if (need_to_add) {
//                     solutions.push_back(cand_sol);
//                 }
//             } else {
//                 // This is not a solution
//                 // We can update the heuristic of this node
//                 vector<size_t> tmp_h = other_node->g + all_pair_lbs[node->id][other_node->id];
//                 NodePtr new_node = make_shared<Node>(node->id, node->g, tmp_h, node->parent, node->path, node->is_forward);

//                 // Check if this lower bound is dominated by the solution
//                 // We cannot remove the node by checking the dominence with nodes in open
//                 bool need_to_add = true;
//                 vector<NodePtr>::iterator it = solutions.begin();
//                 while (it != solutions.end()) {
//                     if (is_dominated_dr(*it, new_node)) {
//                         // one of the solutions is better than the lower bound of the new node
//                         need_to_add = false;
//                         break;
//                     }
//                     ++it;
//                 }

//                 if (need_to_add) {
//                     new_open.push_back(new_node);
//                 }
//             }
//         }
//     }
//     open = new_open;
// }

// void BOPS::update_open(vector<NodePtr>& open, SolutionSet &solutions, size_t target,
//     const vector<NodePtr>& other_open, const vector<NodePtr>& closed, 
//     const vector<NodePtr>& other_closed) {

//     // Debug logs
//     cout << "Start updating open list" << endl;
//     cout << "current open list size: " << open.size() << endl;
//     cout << "current closed list size: " << closed.size() << endl;
//     cout << "other open list size: " << other_open.size() << endl;
//     cout << "other closed list size: " << other_closed.size() << endl;
//     cout << "solution size: " << solutions.size() << endl;

//     // Collect the paths from both open lists
//     if (open.empty()) {
//         cout << "No nodes in open!" << endl;
//         return;
//     } else if (other_open.empty()) {
//         cout << "No nodes in other_open!" << endl;
//         return;
//     }

//     bool is_forward = open[0]->is_forward;

//     list<PathGvalPair> open_paths = get_paths(open);
//     list<PathGvalPair> closed_paths = get_paths(closed);
//     list<PathGvalPair>::iterator op_it = open_paths.begin();
//     while (op_it != open_paths.end()) {
//         if (find(closed_paths.begin(), closed_paths.end(), *op_it) != closed_paths.end()) {
//             op_it = open_paths.erase(op_it);
//         } else {
//             ++ op_it;
//         }
//     }

//     list<PathGvalPair> other_open_paths = get_paths(other_open);
//     list<PathGvalPair> other_closed_paths = get_paths(other_closed);
//     op_it = other_open_paths.begin();
//     while (op_it != other_open_paths.end()) {
//         if (find(other_closed_paths.begin(), other_closed_paths.end(), *op_it) != other_closed_paths.end()) {
//             op_it = other_open_paths.erase(op_it);
//         } else {
//             ++ op_it;
//         }
//     }

//     vector<NodePtr> new_open;
//     for (const auto& p : open_paths) {
//         for (const auto& other_p : other_open_paths) {
//             if (p.first.back() == other_p.first.back()) {
//                 // This is a candidate solution
//                 PathGvalPair _p_ = combine_path_pair(p, other_p, target);
//                 NodePtr can_sol = std::make_shared<Node>(target, _p_.second, 
//                     vector<size_t>(_p_.second.size(), 0));
//                 bool need_to_add = true;
//                 vector<NodePtr>::iterator it = solutions.begin();
//                 while (it != solutions.end()) {
//                     if (is_dominated_dr(can_sol, *it)) {
//                         // candidate node is better than one of the solutions
//                         // TODO: check if this condition ever exists
//                         num_better_sol ++;
//                         it = solutions.erase(it);
//                     } else if (is_dominated_dr(*it, can_sol)) {
//                         // one of the solutions is better than candidate node
//                         num_redundent_sol ++;
//                         need_to_add = false;
//                         break;
//                     } else {
//                         // The cadidate solution is not dominant by the current solution set
//                         ++it;
//                     }
//                 }
//                 if (need_to_add) {
//                     solutions.push_back(can_sol);
//                 }
//             } else {
//                 // This is not a solution
//                 // We can update the heuristic of this node
//                 // We don't need the parent because we already record the path in the node
//                 vector<size_t> tmp_h = other_p.second + all_pair_lbs[p.first.back()][other_p.first.back()];
//                 NodePtr new_node = make_shared<Node>(p.first.back(), p.second, tmp_h, nullptr, is_forward);
//                 new_node->path = p.first;

//                 // Check if this lower bound is dominated by the solution
//                 // We cannot remove the node by checking the dominence with nodes in open
//                 bool need_to_add = true;
//                 vector<NodePtr>::iterator it = solutions.begin();
//                 while (it != solutions.end()) {
//                     if (is_dominated_dr(*it, new_node)) {
//                         // one of the solutions is better than the lower bound of the new node
//                         need_to_add = false;
//                         break;
//                     }
//                     ++it;
//                 }

//                 if (need_to_add) {
//                     new_open.push_back(new_node);
//                 }
//             }
//         }
//     }
//     open = new_open;
// }

void BOPS::update_open(vector<NodePtr>& open, SolutionSet &solutions,  size_t target,
    list<PathGvalPair>& open_paths, const list<PathGvalPair>& other_open_paths) {

    // Debug logs
    // cout << "Start updating open list" << endl;

    bool is_forward = open[0]->is_forward;
    vector<NodePtr> new_open;
    for (const auto& p : open_paths) {
        for (const auto& other_p : other_open_paths) {
            if (p.first.back() == other_p.first.back()) {
                // This is a candidate solution
                PathGvalPair _p_ = combine_path_pair(p, other_p, target);
                NodePtr can_sol = std::make_shared<Node>(target, _p_.second, 
                    vector<size_t>(_p_.second.size(), 0));
                can_sol->path = _p_.first;
                bool need_to_add = true;
                vector<NodePtr>::iterator it = solutions.begin();

                // cout << "_p_ " << can_sol->path.front() << ", " << can_sol->path.back() << endl;
                while (it != solutions.end()) {
                    if (is_bounded(can_sol, *it)) {
                        // one of the solutions is better than candidate node
                        num_redundent_sol ++;
                        need_to_add = false;

                        if (SCREEN) {
                            cout << "worse solution: " << "(" << can_sol->g[0] << ", " << can_sol->g[1]
                            << ") > " << "(" << (*it)->g[0] << " + " << (*it)->h[0] << ", " << 
                            (*it)->g[1] << " + " << (*it)->h[1] << ")" << endl;
                        }

                        break;
                    } else if (is_bounded(*it, can_sol)) {
                        // candidate node is better than one of the solutions
                        if (SCREEN) {
                            cout << "better solution: " << "(" << can_sol->g[0] << ", " << can_sol->g[1]
                            << ") < " << "(" << (*it)->g[0] << " + " << (*it)->h[0] << ", " << 
                            (*it)->g[1] << " + " << (*it)->h[1] << ")" << endl;
                        }

                        num_better_sol ++;
                        it = solutions.erase(it);
                    } else {
                        // The cadidate solution is not dominant by the current solution set
                        ++it;
                    }
                }
                if (need_to_add) {
                    num_sol_middle ++;
                    solutions.push_back(can_sol);
                }
            } else {
                // This is not a solution
                // We can update the heuristic of this node
                // We don't need the parent because we already record the path in the node
                // vector<size_t> tmp_h = other_p.second + all_pair_lbs[p.first.back()][other_p.first.back()];
                vector<size_t> tmp_h = get_diff_heuristic(p.first.back(), other_p.first.back());
                tmp_h = tmp_h + other_p.second;
                NodePtr new_node = make_shared<Node>(p.first.back(), p.second, tmp_h, nullptr, is_forward, &other_p);
                new_node->path = p.first;

                // Check if this lower bound is dominated by the solution
                // We cannot remove the node by checking the dominence with nodes in open
                bool need_to_add = true;
                vector<NodePtr>::iterator it = solutions.begin();
                while (it != solutions.end()) {
                    if (is_bounded(new_node, *it)) {
                        // one of the solutions is better than the lower bound of the new node
                        need_to_add = false;
                        break;
                    }
                    ++it;
                }

                if (need_to_add) {
                    new_open.push_back(new_node);
                }
            }
        }
    }
    open = new_open;
    open_paths.clear();
    cout << "solution size: " << solutions.size() << endl;
}

vector<size_t> BOPS::get_diff_heuristic(size_t loc1, size_t loc2) {
    vector<size_t> h1 = heuristic_f(loc1);
    vector<size_t> h2 = heuristic_f(loc2);
    vector<size_t> tmp_h1 = get_diff(h1, h2);

    h1 = heuristic_b(loc1);
    h2 = heuristic_b(loc2);
    vector<size_t> tmp_h2 = get_diff(h1, h2);
    return get_comax(tmp_h1, tmp_h2);
}