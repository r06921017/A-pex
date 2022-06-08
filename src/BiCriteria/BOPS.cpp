#include <memory>
#include <algorithm>
#include <time.h>

#include "BOPS.h"

#define SCREEN true

BOPS::BOPS(const AdjacencyMatrix &adj_matrix, Pair<double> eps, Heuristic &h_f, Heuristic &h_b,
    const LoggerPtr logger, size_t lh_f, size_t lh_b) : 
    BOAStar(adj_matrix, eps, logger), lookahead_f(lh_f), lookahead_b(lh_b),
    heuristic_f(h_f), heuristic_b(h_b) {}

// TODO: only count the 1st element of the front to front heuristic, and count the 2nd if needed
void BOPS::operator() (size_t source, size_t target, SolutionSet & solutions,
    unsigned int time_limit) {

    start_time = std::clock();
    this->start_logging(source, target);

    NodePtr node_f;
    NodePtr node_b;
    NodePtr next;

    // Saving all the unused NodePtrs in a vector improves performace for some reason
    vector<NodePtr> closed_f;  // forward
    vector<NodePtr> closed_b;  // backward
    vector<NodePtr> candidate_sols;  // set for the condidate solutions, need to concadinate with either the forward or backward open list

    // Vector to hold mininum cost of 2nd criteria per node
    vector<size_t> min_g2_f(this->adj_matrix.size()+1, MAX_COST);
    vector<size_t> min_g2_b(this->adj_matrix.size()+1, MAX_COST);

    // Init open heap
    Node::more_than_full_cost more_than;
    vector<NodePtr> open_f;
    make_heap(open_f.begin(), open_f.end(), more_than);
    vector<NodePtr> open_b;
    make_heap(open_b.begin(), open_b.end(), more_than);

    node_f = make_shared<Node>(source, std::vector<size_t>(2,0), heuristic_f(source), nullptr, true);
    node_f->set_path();

    node_b = make_shared<Node>(target, std::vector<size_t>(2,0), heuristic_b(target), nullptr, false);
    node_b->set_path();

    node_f->set_h_node(node_b);
    node_b->set_h_node(node_f);

    open_f.push_back(node_f);
    std::push_heap(open_f.begin(), open_f.end(), more_than);
    list<PathGvalPair> open_f_paths = get_paths(open_f);

    open_b.push_back(node_b);
    std::push_heap(open_b.begin(), open_b.end(), more_than);
    list<PathGvalPair> open_b_paths = get_paths(open_b);

    while (!open_f.empty()) {
        // Begin with forward search
        if (SCREEN)
            cout << "Start forward search with open_f size: " << open_f.size() << endl;

        // Update the heuristic for the backward open list given the forward open list
        for (auto& _node_ : open_f) {
            list<HeuristicNodePair> tmp_h_list;
            for (const auto& _other_node_ : open_b) {
                if (_other_node_->is_cand) continue;
                vector<size_t> tmp_h = get_diff_heuristic(_node_->id, _other_node_->id);  // f2f
                tmp_h += _other_node_->g;
                bool is_b = false;
                list<HeuristicNodePair>::iterator h_it = tmp_h_list.begin();
                while (h_it != tmp_h_list.end()) {
                    if (is_bounded(h_it->first, tmp_h)) {
                        h_it = tmp_h_list.erase(h_it);
                    } else if (is_bounded(tmp_h, h_it->first)) {
                        is_b = true;
                        break;
                    } else {
                        ++ h_it;
                    }
                }

                if (!is_b) {
                    tmp_h_list.emplace_back(tmp_h, _other_node_);
                    // if (SCREEN) {
                    //     vector<size_t> h1_f = heuristic_f(_node_->id);
                    //     vector<size_t> h1_b = heuristic_b(_node_->id);
                    //     cout << "curr node id: " << _node_->id << endl;
                    //     cout << "\th_f: [" << h1_f[0] << ", " << h1_f[1] << "], h_b: [" << h1_b[0] 
                    //         << ", " << h1_b[1] << "]" << endl;
                    //     h1_f = heuristic_f(_other_node_->id);
                    //     h1_b = heuristic_b(_other_node_->id);
                    //     cout << "other node id: " << _other_node_->id << endl;
                    //     cout << "\th_f: [" << h1_f[0] << ", " << h1_f[1] << "], h_b: [" << h1_b[0] 
                    //         << ", " << h1_b[1] << "]" << endl;
                    //     cout << "f2f h: [" << tmp_h[0] << ", " << tmp_h[1] << "]" << endl;
                    //     cout << endl;
                    // }
                }
            }
            _node_->set_h_val(tmp_h_list);

            // if (SCREEN) {
            //     _node_->print_h();
            //     cout << endl;
            // }
        }

        if (!candidate_sols.empty()) {
            if (SCREEN) {
                cout << "++++++++++++++++++++" << endl;
                cout << "Before adding candidate solutions:" << endl;
                print_list(open_f, more_than);
                cout << "open_f.size() = " << open_f.size() << endl;
            }

            open_f.insert(open_f.end(), candidate_sols.begin(), candidate_sols.end());
            make_heap(open_f.begin(), open_f.end(), more_than);

            if (SCREEN) {
                cout << "After adding candidate solutions: " << endl;
                print_list(open_f, more_than);
                cout << "open_f.size() = " << open_f.size() << endl;
                cout << "candidate_sols.size() = " << candidate_sols.size() << endl;
                cout << "++++++++++++++++++++" << endl;
            }
        }

        while (!open_f.empty()) {
            if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit) {
                this->end_logging(solutions, false);
                return;
            } else if (num_expansion_f > lookahead_f) {  // Switch to backward search
                num_expansion += num_expansion_f;
                num_generation += num_generation_f;
                num_expansion_f = 0;
                num_generation_f = 0;

                if (SCREEN) {
                    pop_heap(open_f.begin(), open_f.end(), more_than);
                    NodePtr tmp_node_f = open_f.back();
                    cout << "\nForward search end" << endl;
                    cout << "Current open_f top: " << *tmp_node_f << endl;
                    cout << "open_f size: " << open_f.size() << endl;
                    cout << "solutions size: " << solutions.size() << endl;
                    cout << "candidate_sols size: " << candidate_sols.size() << endl;
                }
                break;
            }

            // Pop min from queue and process
            pop_heap(open_f.begin(), open_f.end(), more_than);
            node_f = open_f.back();
            open_f.pop_back();
            // assert(node_f->path.back() == node_f->id);

            // Dominance check
            // if ((!node_f->is_cand && ((1+this->eps[1])*node_f->f[1]) >= min_g2_f[target]) ||
            //     (!node_f->is_cand && node_f->g[1] >= min_g2_f[node_f->id]) || 
            //     (node_f->is_cand && node_f->g[1] >= min_g2_f[target])) {
            //     reinsert(node_f, open_f, closed_f, more_than);
            //     continue;
            // }
            if ((((1+this->eps[1])*node_f->f[1]) >= min_g2_f[target]) ||
                (node_f->g[1] >= min_g2_f[node_f->id])) {
                if (node_f->is_cand) {
                    cout << "can_sol is dominant" << endl;
                    cout << *node_f << endl;
                }
                reinsert(node_f, open_f, closed_f, more_than);
                continue;
            }

            // We decide to expand this node
            // TODO: might need the other heuristics for the pathmax heuristic 
            // if front to front is too time-consuming
            if (SCREEN) {
                cout << "---------------------------------" << endl;
                cout << "expand node_f" << endl;
                cout << *node_f << endl;
                cout << "---------------------------------" << endl;
            }

            min_g2_f[node_f->id] = node_f->g[1];
            num_expansion_f += 1;

            // Find one solution during the search
            if (node_f->is_cand) {
                solutions.push_back(node_f);
                num_sol_front ++;
                log_solution(node_f);

                if (SCREEN) {
                    cout << "This is a solution!" << endl;
                    cout << *node_f;
                    cout << endl;
                    cout << "----- open_f -----" << endl;
                    print_list(open_f, more_than, 3);
                }

                // This is a solution found in the previous iterations.
                // Remove it from the candidate solution set
                assert(node_f->h == (size_t) 0);
                    // cout << "candidate_sols: " << candidate_sols.size() << endl;
                candidate_sols.erase(remove(candidate_sols.begin(), candidate_sols.end(), node_f), 
                    candidate_sols.end());
                    // cout << "new candidate_sols: " << candidate_sols.size() << endl;

                // Push to the closed list
                min_g2_f[target] = node_f->g[1];  // update for global dominance check
                closed_f.push_back(node_f);
                continue;
            } else if (node_f->id == node_f->h_node->id) {
                // This is a candidate solution, will be checked in the opposite search
                // Don't need to put in the node to the OPEN
                NodePtr can_sol = make_shared<Node>(node_f->id, node_f->g + node_f->h_node->g, 
                    vector<size_t>(node_f->g.size(), 0), node_f->parent, true, node_f->h_node, true);
                can_sol->combine_path();
                candidate_sols.push_back(can_sol);
                open_f.push_back(can_sol);
                push_heap(open_f.begin(), open_f.end(), more_than);

                if (SCREEN) {
                    cout << "This is a can_sol" << endl;
                    cout << *can_sol;
                    cout << endl;
                    cout << "----- open_f -----" << endl;
                    print_list(open_f, more_than, 3);
                }

                // Reinsert the node with another heuristic to the open
                reinsert(node_f, open_f, closed_f, more_than);
                continue;
            }

            // Check to which neighbors we should extend the paths
            const vector<Edge> &outgoing_edges = adj_matrix[node_f->id];
            if (SCREEN)
                cout << "outgoing_edges size: " << outgoing_edges.size() << endl;
            for (auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                if (node_f->parent != nullptr && p_edge->target == node_f->parent->id) 
                    continue;
                size_t next_id = p_edge->target;
                vector<size_t> next_g = node_f->g + p_edge->cost;

                list<HeuristicNodePair> h_list;
                for (const auto& other_node : open_b) {
                    vector<size_t> next_h = get_diff_heuristic(next_id, other_node->id);  // f2f
                    next_h += other_node->g;

                    // Dominance check
                    if ((((1+this->eps[1])*(next_g[1]+next_h[1])) >= min_g2_f[target]) ||
                        (next_g[1] >= min_g2_f[next_id])) {
                        // if (SCREEN)
                        //     cout << "is dominated in " << next_id << endl;    
                        continue;
                    }

                    // Check if the heuristic is dominant
                    bool need_continue = false;
                    list<HeuristicNodePair>::iterator h_it = h_list.begin();
                    while (h_it != h_list.end()) {
                        if (is_bounded(h_it->first, next_h)) {
                            h_it = h_list.erase(h_it);
                        } else if (is_bounded(next_h, h_it->first)) {
                            need_continue = true;
                            break;
                        } else {
                            ++ h_it;
                        }
                    }
                    if (need_continue) continue;
                    h_list.emplace_back(next_h, other_node);
                }

                if (!h_list.empty()) {
                    // If not dominated create node and push to queue
                    // Creation is defered after dominance check as it is
                    // relatively computational heavy and should be avoided if possible
                    next = make_shared<Node>(next_id, next_g, h_list, node_f, node_f->is_forward);
                    next->set_path();
                    open_f.push_back(next);
                    push_heap(open_f.begin(), open_f.end(), more_than);
                    num_generation_f +=1;
                }

                if (SCREEN) {
                    cout << "\tgenerate node" << endl;
                    cout << "\t" << *next << endl;
                }
            }
            closed_f.push_back(node_f);  // We fully expand this node, so add to closed list
            if (SCREEN) {
                cout << "generation done, open_f.size: " << open_f.size() << ", closed_f.size: "
                    << closed_f.size() << endl;
                cout << endl;
            }
        }

        // Backward search
        if (SCREEN)
            cout << "Start backward search with open_b size: " << open_b.size() << endl;

        // Update the heuristic for the backward open list given the forward open list
        for (auto& _node_ : open_b) {
            list<HeuristicNodePair> tmp_h_list;
            for (const auto& _other_node_ : open_f) {
                if (_other_node_->is_cand) continue;
                vector<size_t> tmp_h = get_diff_heuristic(_node_->id, _other_node_->id);  // f2f
                tmp_h += _other_node_->g;
                bool is_b = false;
                list<HeuristicNodePair>::iterator h_it = tmp_h_list.begin();
                while (h_it != tmp_h_list.end()) {
                    if (is_bounded(h_it->first, tmp_h)) {
                        h_it = tmp_h_list.erase(h_it);
                    } else if (is_bounded(tmp_h, h_it->first)) {
                        is_b = true;
                        break;
                    } else {
                        ++ h_it;
                    }
                }

                if (!is_b) {
                    tmp_h_list.emplace_back(tmp_h, _other_node_);
                    // if (SCREEN) {
                    //     vector<size_t> h1_f = heuristic_f(_node_->id);
                    //     vector<size_t> h1_b = heuristic_b(_node_->id);
                    //     cout << "curr node id: " << _node_->id << endl;
                    //     cout << "\th_f: [" << h1_f[0] << ", " << h1_f[1] << "], h_b: [" << h1_b[0] 
                    //         << ", " << h1_b[1] << "]" << endl;
                    //     h1_f = heuristic_f(_other_node_->id);
                    //     h1_b = heuristic_b(_other_node_->id);
                    //     cout << "other node id: " << _other_node_->id << endl;
                    //     cout << "\th_f: [" << h1_f[0] << ", " << h1_f[1] << "], h_b: [" << h1_b[0] 
                    //         << ", " << h1_b[1] << "]" << endl;
                    //     cout << "f2f h: [" << tmp_h[0] << ", " << tmp_h[1] << "]" << endl;
                    //     cout << endl;
                    // }
                }
            }
            _node_->set_h_val(tmp_h_list);

            if (SCREEN) {
                cout << *_node_ << endl;
                _node_->print_h();
                // cout << _node_->other_h.size();
                cout << endl;
            }
        }

        // if (!candidate_sols.empty()) {
        //     if (SCREEN) {
        //         cout << "++++++++++++++++++++" << endl;
        //         cout << "Before adding candidate solutions:" << endl;
        //         print_list(open_b, more_than);
        //         cout << "open_b.size() = " << open_b.size() << endl;
        //     }

        //     // for (const auto& can_sol : candidate_sols) {
        //     //     open_b.push_back(can_sol);
        //     //     push_heap(open_b.begin(), open_b.end(), more_than);
        //     // }
        //     open_b.insert(open_b.end(), candidate_sols.begin(), candidate_sols.end());
        //     if (SCREEN) {
        //         cout << "open_b before sorting ..." << endl;
        //         print_list(open_b, more_than);
        //         cout << endl;
        //     }
        //     make_heap(open_b.begin(), open_b.end(), more_than);

        //     if (SCREEN) {
        //         cout << "After adding candidate solutions: " << endl;
        //         print_list(open_b, more_than);
        //         cout << "open_b.size() = " << open_b.size() << endl;
        //         cout << "candidate_sols.size() = " << candidate_sols.size() << endl;
        //         cout << "++++++++++++++++++++" << endl;
        //     }
        // }

        while (!open_b.empty()) {
            if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit) {
                this->end_logging(solutions, false);
                return;
            } else if (num_expansion_b > lookahead_b) {  // Switch to forward search
                num_expansion += num_expansion_b;
                num_generation += num_generation_b;
                num_expansion_b = 0;
                num_generation_b = 0;

                if (SCREEN) {
                    pop_heap(open_b.begin(), open_b.end(), more_than);
                    NodePtr tmp_node_b = open_b.back();
                    cout << "\nBackward search end" << endl;
                    cout << "Current open_b top: " << *tmp_node_b << endl;
                    cout << "open_b size: " << open_b.size() << endl;
                    cout << "solutions size: " << solutions.size() << endl;
                    cout << "candidate_sols size: " << candidate_sols.size() << endl;
                }
                break;
            }

            // Pop min from queue and process
            std::pop_heap(open_b.begin(), open_b.end(), more_than);
            node_b = open_b.back();
            open_b.pop_back();
            assert(!node_b->is_cand && !node_b->is_forward);
            // assert(node_b->path.back() == node_b->id);

            if (SCREEN) {
                cout << "---------------------------------" << endl;
                cout << "expand node_b" << endl;
                cout << *node_b << endl;
                cout << "min_g2_b[source]: " << min_g2_b[source] << endl;
                cout << "min_g2_b[node_b->id]: " << min_g2_b[node_b->id] << endl;
                cout << "---------------------------------" << endl;
            }

            // Dominance check
            if ((((1+this->eps[1])*node_b->f[1]) >= min_g2_b[source]) ||
                (node_b->g[1] >= min_g2_b[node_b->id])) {
                reinsert(node_b, open_b, closed_b, more_than);
                continue;
            }

            // We decide to expand this node
            // TODO: might need the other heuristics for the pathmax heuristic 
            // if front to front is too time-consuming
            min_g2_b[node_b->id] = node_b->g[1];
            num_expansion_b += 1;

            // Find one solution during the search
            if (node_b->id == node_b->h_node->id) {
                // This is a candidate solution, will be checked in the opposite search
                // Don't need to put in the node to the OPEN
                NodePtr can_sol = make_shared<Node>(node_b->id, node_b->g + node_b->h_node->g, 
                    vector<size_t>(node_b->g.size(), 0), node_b->parent, false, node_b->h_node, true);
                can_sol->combine_path();
                candidate_sols.push_back(can_sol);

                if (SCREEN) {
                    cout << "This is a can_sol" << endl;
                    cout << *can_sol;
                    cout << endl;
                }

                // Reinsert the node with another heuristic to the open
                reinsert(node_b, open_b, closed_b, more_than);
                if (SCREEN) {
                    cout << "node_b" << endl;
                    cout << *node_b << endl;
                    cout << endl;
                }
                continue;
            }

            // Check to which neighbors we should extend the paths
            const vector<Edge> &outgoing_edges = adj_matrix[node_b->id];
            if (SCREEN)
                cout << "outgoing_edges size: " << outgoing_edges.size() << endl;
            for (auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                if (node_b->parent != nullptr && p_edge->target == node_b->parent->id) 
                    continue;
                size_t next_id = p_edge->target;
                vector<size_t> next_g = node_b->g + p_edge->cost;

                list<HeuristicNodePair> h_list;
                for (const auto& other_node : open_f) {
                    vector<size_t> next_h = get_diff_heuristic(next_id, other_node->id);  // f2f
                    next_h += other_node->g;

                    // Dominance check
                    if ((((1+this->eps[1])*(next_g[1]+next_h[1])) >= min_g2_b[source]) ||
                        (next_g[1] >= min_g2_b[next_id])) {
                        // if (SCREEN)
                        //     cout << "is dominated in " << next_id << endl;    
                        continue;
                    }

                    // Check if the heuristic is dominant
                    bool need_continue = false;
                    list<HeuristicNodePair>::iterator h_it = h_list.begin();
                    while (h_it != h_list.end()) {
                        if (is_bounded(h_it->first, next_h)) {
                            h_it = h_list.erase(h_it);
                        } else if (is_bounded(next_h, h_it->first)) {
                            need_continue = true;
                            break;
                        } else {
                            ++ h_it;
                        }
                    }
                    if (need_continue) continue;
                    h_list.emplace_back(next_h, other_node);
                }

                if (!h_list.empty()) {
                    // If not dominated create node and push to queue
                    // Creation is defered after dominance check as it is
                    // relatively computational heavy and should be avoided if possible
                    next = make_shared<Node>(next_id, next_g, h_list, node_b, node_b->is_forward);
                    next->set_path();
                    open_b.push_back(next);
                    push_heap(open_b.begin(), open_b.end(), more_than);
                    num_generation_b +=1;
                }

                if (SCREEN) {
                    cout << "\tgenerate node" << endl;
                    cout << "\t" << *next << endl;
                    cout << "\tnumber of h: " << next->other_h.size() + 1 << endl; 
                    // next->print_h();
                }
            }
            closed_b.push_back(node_b);  // We fully expand this node, so add to closed list
            if (SCREEN) {
                cout << "generation done, open_b.size: " << open_b.size() << ", closed_b.size: "
                    << closed_b.size() << endl;
                cout << endl;
            }
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

// void BOPS::update_open(vector<NodePtr>& open, SolutionSet &solutions,  size_t target,
//     list<PathGvalPair>& open_paths, const list<PathGvalPair>& other_open_paths) {

//     // Debug logs
//     // cout << "Start updating open list" << endl;

//     bool is_forward = open[0]->is_forward;
//     vector<NodePtr> new_open;
//     for (const auto& p : open_paths) {
//         for (const auto& other_p : other_open_paths) {
//             if (p.first.back() == other_p.first.back()) {
//                 // This is a candidate solution
//                 PathGvalPair _p_ = combine_path_pair(p, other_p, target);
//                 NodePtr can_sol = std::make_shared<Node>(target, _p_.second, 
//                     vector<size_t>(_p_.second.size(), 0));
//                 can_sol->path = _p_.first;
//                 bool need_to_add = true;
//                 vector<NodePtr>::iterator it = solutions.begin();

//                 // cout << "_p_ " << can_sol->path.front() << ", " << can_sol->path.back() << endl;
//                 while (it != solutions.end()) {
//                     if (is_bounded(can_sol, *it)) {
//                         // one of the solutions is better than candidate node
//                         num_redundent_sol ++;
//                         need_to_add = false;

//                         if (SCREEN) {
//                             cout << "worse solution: " << "(" << can_sol->g[0] << ", " << can_sol->g[1]
//                             << ") > " << "(" << (*it)->g[0] << " + " << (*it)->h[0] << ", " << 
//                             (*it)->g[1] << " + " << (*it)->h[1] << ")" << endl;
//                         }

//                         break;
//                     } else if (is_bounded(*it, can_sol)) {
//                         // candidate node is better than one of the solutions
//                         if (SCREEN) {
//                             cout << "better solution: " << "(" << can_sol->g[0] << ", " << can_sol->g[1]
//                             << ") < " << "(" << (*it)->g[0] << " + " << (*it)->h[0] << ", " << 
//                             (*it)->g[1] << " + " << (*it)->h[1] << ")" << endl;
//                         }

//                         num_better_sol ++;
//                         it = solutions.erase(it);
//                     } else {
//                         // The cadidate solution is not dominant by the current solution set
//                         ++it;
//                     }
//                 }
//                 if (need_to_add) {
//                     num_sol_middle ++;
//                     solutions.push_back(can_sol);
//                 }
//             } else {
//                 // This is not a solution
//                 // We can update the heuristic of this node
//                 // We don't need the parent because we already record the path in the node
//                 // vector<size_t> tmp_h = other_p.second + all_pair_lbs[p.first.back()][other_p.first.back()];
//                 vector<size_t> tmp_h = get_diff_heuristic(p.first.back(), other_p.first.back());
//                 tmp_h = tmp_h + other_p.second;
//                 NodePtr new_node = make_shared<Node>(p.first.back(), p.second, tmp_h, nullptr, is_forward, &other_p);
//                 new_node->path = p.first;

//                 // Check if this lower bound is dominated by the solution
//                 // We cannot remove the node by checking the dominence with nodes in open
//                 bool need_to_add = true;
//                 vector<NodePtr>::iterator it = solutions.begin();
//                 while (it != solutions.end()) {
//                     if (is_bounded(new_node, *it)) {
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
//     open_paths.clear();
//     cout << "solution size: " << solutions.size() << endl;
// }

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
        cout << "push to open" << endl; 
    } else {
        closed.push_back(node);
        cout << "push to closed" << endl; 
    }
}
