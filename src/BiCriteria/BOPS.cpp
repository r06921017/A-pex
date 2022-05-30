#include <memory>
#include <algorithm>
#include <time.h>

#include "BOPS.h"

BOPS::BOPS(const AdjacencyMatrix &adj_matrix, Pair<double> eps, const LoggerPtr logger, size_t lh_f,
    size_t lh_b) : BOAStar(adj_matrix, eps, logger), lookahead_f(lh_f), lookahead_b(lh_b) {}

void BOPS::operator() (size_t source, size_t target, Heuristic &heuristic, 
    Heuristic &heuristic_b, SolutionSet & solutions, AdjacencyMatrix& graph, size_t graph_size, 
    unsigned int time_limit) {

    // Get the all pair shortest path heuristic, which is a set of vectors
    all_pair_lbs = std::vector<std::vector<std::vector<size_t>>>(graph_size);
    for (size_t i = 0; i < graph_size; i++) {
        all_pair_lbs[i] = std::vector<std::vector<size_t>>(graph_size);
    }

    for (size_t loc1=0; loc1 < graph_size; loc1++) {
        ShortestPathHeuristic tmp_sp_heuristic(loc1, graph_size, graph);
        for (size_t loc2=0; loc2 < graph_size; loc2++) {
            all_pair_lbs[loc1][loc2] = tmp_sp_heuristic(loc2);
        }
    }

    start_time = std::clock();
    this->start_logging(source, target);

    NodePtr node_f;
    NodePtr node_b;
    NodePtr next;
    list<NodePtr> tmp_path;

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

    node_f = std::make_shared<Node>(source, std::vector<size_t>(2,0), heuristic(source));
    node_f->set_path();
    open_f.push_back(node_f);
    std::push_heap(open_f.begin(), open_f.end(), more_than);

    node_b = std::make_shared<Node>(target, std::vector<size_t>(2,0), heuristic_b(target));
    node_f->set_path();
    open_b.push_back(node_b);
    std::push_heap(open_b.begin(), open_b.end(), more_than);

    while (!open_f.empty() && !open_b.empty())
    {
        // Begin with forward search
        update_open(open_f, open_b, target, solutions);
        while (!open_f.empty()) {
            if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit){
                this->end_logging(solutions, false);
                return;
            }

            else if (num_expansion_f > lookahead_f)  // Switch to backward search
            {
                num_expansion += num_expansion_f;
                num_generation += num_generation_f;
                num_expansion_f = 0;
                num_generation_f = 0;
                break;
            }

            // Pop min from queue and process
            std::pop_heap(open_f.begin(), open_f.end(), more_than);
            node_f = open_f.back();
            open_f.pop_back();
            num_generation_f +=1;

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
                solutions.push_back(node_f);
                log_solution(node_f);
                continue;
            }

            // Check to which neighbors we should extend the paths
            const std::vector<Edge> &outgoing_edges = adj_matrix[node_f->id];
            for (auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                size_t next_id = p_edge->target;
                std::vector<size_t> next_g = {node_f->g[0]+p_edge->cost[0], 
                    node_f->g[1]+p_edge->cost[1]};
                std::vector<size_t> next_h = heuristic(next_id);

                // Dominance check
                if ((((1+this->eps[1])*(next_g[1]+next_h[1])) >= min_g2_f[target]) ||
                    (next_g[1] >= min_g2_f[next_id])) {
                    continue;
                }

                // If not dominated create node and push to queue
                // Creation is defered after dominance check as it is
                // relatively computational heavy and should be avoided if possible
                next = std::make_shared<Node>(next_id, next_g, next_h, node_f);
                next->set_path();

                open_f.push_back(next);
                std::push_heap(open_f.begin(), open_f.end(), more_than);

                closed_f.push_back(node_f);
            }
        }

        // Backward search
        update_open(open_b, open_f, target, solutions);
        while (!open_b.empty()) {
            if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit){
                this->end_logging(solutions, false);
                return;
            }

            else if (num_expansion_b > lookahead_b)  // Switch to forward search
            {
                num_expansion += num_expansion_b;
                num_generation += num_generation_b;
                num_expansion_b = 0;
                num_generation_b = 0;
                break;
            }

            // Pop min from queue and process
            std::pop_heap(open_b.begin(), open_b.end(), more_than);
            node_b = open_b.back();
            open_b.pop_back();
            num_generation_b +=1;
            assert(node_b->path.back() == node_b->id);

            // Dominance check
            if ((((1+this->eps[1])*node_b->f[1]) >= min_g2_b[source]) ||
                (node_b->g[1] >= min_g2_b[node_b->id])) {
                closed_b.push_back(node_b);
                continue;
            }

            min_g2_b[node_b->id] = node_b->g[1];
            num_expansion += 1;

            // Find one solution during the search
            // We already did this in update_open function
            // if (node_b->id == source) {
            //     solutions.push_back(node_b);
            //     log_solution(node_b);
            //     continue;
            // }

            // Check to which neighbors we should extend the paths
            const std::vector<Edge> &outgoing_edges = adj_matrix[node_f->id];
            for (auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                size_t next_id = p_edge->target;
                std::vector<size_t> next_g = {node_f->g[0]+p_edge->cost[0], 
                    node_f->g[1]+p_edge->cost[1]};
                std::vector<size_t> next_h = heuristic(next_id);

                // Dominance check
                if ((((1+this->eps[1])*(next_g[1]+next_h[1])) >= min_g2_f[target]) ||
                    (next_g[1] >= min_g2_f[next_id])) {
                    continue;
                }

                // If not dominated create node and push to queue
                // Creation is defered after dominance check as it is
                // relatively computational heavy and should be avoided if possible
                next = std::make_shared<Node>(next_id, next_g, next_h, node_b, node_b->is_forward);
                next->set_path();
                open_b.push_back(next);
                std::push_heap(open_b.begin(), open_b.end(), more_than);
                closed_f.push_back(node_f);
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

void BOPS::update_open(vector<NodePtr>& open, const vector<NodePtr>& other_open,
    size_t target, SolutionSet &solutions) {

    // Collect the paths from both open lists
    if (open.empty()) {
        cout << "No nodes in open!" << endl;
        return;
    } else if (other_open.empty()) {
        cout << "No nodes in other_open!" << endl;
        return;
    }

    bool is_forward = open[0]->is_forward;
    list<PathGvalPair> open_paths = get_paths(open);
    list<PathGvalPair> other_open_paths = get_paths(other_open);
    vector<NodePtr> new_open;

    for (const auto& p : open_paths) {
        for (const auto& other_p : other_open_paths) {
            if (p.first.back() == other_p.first.back()) {
                // This is a candidate solution
                PathGvalPair _p_ = combine_path_pair(p, other_p, target);
                NodePtr can_sol = std::make_shared<Node>(target, _p_.second, 
                    vector<size_t>(_p_.second.size(), 0));
                bool need_to_add = true;
                vector<NodePtr>::iterator it = solutions.begin();
                while (it != solutions.end()) {
                    if (is_dominated_dr(can_sol, *it)) {
                        // candidate node is better than one of the solutions
                        // TODO: check if this condition ever exists
                        num_better_sol ++;
                        it = solutions.erase(it);
                    } else if (is_dominated_dr(*it, can_sol)) {
                        // one of the solutions is better than candidate node
                        num_redundent_sol ++;
                        need_to_add = false;
                        break;
                    } else {
                        // The cadidate solution is not dominant by the current solution set
                        ++it;
                    }
                }
                if (need_to_add) {
                    solutions.push_back(can_sol);
                }
            } else {
                // This is not a solution
                // We can update the heuristic of this node
                // We don't need the parent because we already record the path in the node
                vector<size_t> tmp_h = other_p.second + all_pair_lbs[p.first.back()][other_p.first.back()];
                NodePtr new_node = make_shared<Node>(p.first.back(), p.second, tmp_h, nullptr, is_forward);

                // Check if this lower bound is dominated by the solution
                // We cannot remove the node by checking the dominence with nodes in open
                bool need_to_add = true;
                vector<NodePtr>::iterator it = solutions.begin();
                while (it != solutions.end()) {
                    if (is_dominated_dr(*it, new_node)) {
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
}