#include <memory>
#include <algorithm>
#include <time.h>
#include <iomanip>

#include "BOAStar.h"

BOAStar::BOAStar(const AdjacencyMatrix &adj_matrix, Pair<double> eps, 
    const LoggerPtr logger, uint screen) :
    AbstractSolver(adj_matrix, {eps[0], eps[1]}, logger, screen) {}

void BOAStar::operator()(size_t source, size_t target, Heuristic &heuristic, 
    SolutionSet &solutions, uint time_limit) {
    clock_t tmp_start;
    start_time = clock();
    this->start_logging(source, target);

    NodePtr node;
    NodePtr next;

    // Saving all the unused NodePtrs in a vector improves performace for some reason
    vector<NodePtr> closed;

    // Vector to hold mininum cost of 2nd criteria per node
    vector<size_t> min_g2(this->adj_matrix.size()+1, MAX_COST);

    // Init open heap
    Node::more_than_full_cost more_than;
    vector<NodePtr> open;
    make_heap(open.begin(), open.end(), more_than);

    node = make_shared<Node>(source, vector<size_t>(2,0), heuristic(source));
    open.push_back(node);
    push_heap(open.begin(), open.end(), more_than);

    while (open.empty() == false) {
        if ((clock() - start_time)/CLOCKS_PER_SEC > time_limit){
            if (screen > DEBUG_LOG_EXPANSION) this->end_logging(solutions, false);
            return;
        }

        // Pop min from queue and process
        pop_heap(open.begin(), open.end(), more_than);
        node = open.back();
        open.pop_back();

        // Dominance check
        tmp_start = clock();
        if ((((1+this->eps[1])*node->f[1]) >= min_g2[target]) ||
            (node->g[1] >= min_g2[node->id])) {
            closed.push_back(node);
            runtime_global_check += clock() - tmp_start;
            continue;
        }
        runtime_global_check += clock() - tmp_start;

        min_g2[node->id] = node->g[1];
        num_expansion += 1;

        if (node->id == target) {
            solutions.push_back(node);
            log_solution(node);

            // cout << "g value: " << node->g[0] << ", " << node->g[1] << endl;
            // cout << "f value: " << node->f[0] << ", " << node->f[1] << endl;
            // NodePtr tmp_n = node;
            // while (tmp_n->parent != nullptr)
            // {
            //     cout << right << setw(5) << tmp_n->id << ",";
            //     tmp_n = tmp_n->parent;
            // }
            // cout << right << setw(5) << source << endl;

            // tmp_n = node;
            // while (tmp_n->parent != nullptr)
            // {
            //     cout << right << setw(5) << tmp_n->g[0] << ",";
            //     tmp_n = tmp_n->parent;
            // }
            // cout << right << setw(5) << 0 << endl;

            // tmp_n = node;
            // while (tmp_n->parent != nullptr)
            // {
            //     cout << right << setw(5) << tmp_n->g[1] << ",";
            //     tmp_n = tmp_n->parent;
            // }
            // cout << right << setw(5) << 0 << endl;

            continue;
        }

        // Check to which neighbors we should extend the paths
        const vector<Edge> &outgoing_edges = adj_matrix[node->id];
        branching_factor_f += outgoing_edges.size() - 1;
        if (screen > DEBUG_LOG_DETAILES)
            cout << "\toutgoing_edges size: " << outgoing_edges.size() << endl;
        for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            // Ignore the node we are coming from
            if (node->parent != nullptr && p_edge->target == node->parent->id) continue;
            clock_t node_gen_start = clock();
            size_t next_id = p_edge->target;
            vector<size_t> next_g = {node->g[0]+p_edge->cost[0], node->g[1]+p_edge->cost[1]};
            auto next_h = heuristic(next_id);

            // Dominance check
            tmp_start = clock();
            if ((((1+this->eps[1])*(next_g[1]+next_h[1])) >= min_g2[target]) ||
                (next_g[1] >= min_g2[next_id])) {
                runtime_global_check += clock() - tmp_start;
                continue;
            }
            runtime_global_check += clock() - tmp_start;

            // If not dominated create node and push to queue
            // Creation is defered after dominance check as it is
            // relatively computational heavy and should be avoided if possible
            next = make_shared<Node>(next_id, next_g, next_h, node);
            num_generation +=1;

            open.push_back(next);
            push_heap(open.begin(), open.end(), more_than);
            runtime_node_gen += clock() - node_gen_start;
        }
        closed.push_back(node);
    }

    if (screen > DEBUG_LOG_EXPANSION) this->end_logging(solutions);
}


inline bool is_dominated(SolutionSet & solutions, NodePtr node, Pair<double> eps = {0,0}){
    for (auto sol: solutions){
        if (sol->g[0] <= (1 + eps[0]) * node->g[0] + node->h[0] &&
            sol->g[1] <= (1 + eps[1]) * node->g[1] + node->h[1]
            ){
            return true;
        }
    }
    return false;
}




void BOAStar::log_solution(NodePtr node){
    solution_log.push_back({clock() - start_time, node});
}
