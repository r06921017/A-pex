#include <memory>
#include <algorithm>
#include <time.h>

#include "BOZStar.h"

BOZStar::BOZStar(const AdjacencyMatrix &adj_matrix, Pair<double> eps, const LoggerPtr logger, int lh_f, int lh_b) :
    BOAStar(adj_matrix, eps, logger), lookahead_f(lh_f), lookahead_b(lh_b) {}

void BOZStar::operator() (size_t source, size_t target, Heuristic &heuristic, 
    Heuristic &heuristic_b, SolutionSet & solutions, AdjacencyMatrix& graph, size_t graph_size, 
    unsigned int time_limit) {

    // Get the all pair shortest path heuristic, which is a set of vectors
    all_pair_path_lbs = std::vector<std::vector<std::vector<size_t>>>(graph_size);
    for (size_t i = 0; i < graph_size; i++) {
        all_pair_path_lbs[i] = std::vector<std::vector<size_t>>(graph_size);
    }

    for (size_t loc1=0; loc1 < graph_size; loc1++) {
        ShortestPathHeuristic tmp_sp_heuristic(loc1, graph_size, graph);
        for (size_t loc2=0; loc2 < graph_size; loc2++) {
            all_pair_path_lbs[loc1][loc2] = tmp_sp_heuristic(loc2);
        }
    }

    start_time = std::clock();
    this->start_logging(source, target);

    NodePtr node_f;
    NodePtr node_b;
    NodePtr next;

    // Saving all the unused NodePtrs in a vector improves performace for some reason
    std::vector<NodePtr> closed_f;  // forward
    std::vector<NodePtr> closed_b;  // backward

    // Vector to hold mininum cost of 2nd criteria per node
    std::vector<size_t> min_g2_f(this->adj_matrix.size()+1, MAX_COST);
    std::vector<size_t> min_g2_b(this->adj_matrix.size()+1, MAX_COST);

    // Init open heap
    Node::more_than_full_cost more_than_f;
    std::vector<NodePtr> open_f;
    std::make_heap(open_f.begin(), open_f.end(), more_than_f);

    Node::more_than_full_cost_b more_than_b;
    std::vector<NodePtr> open_b;
    std::make_heap(open_b.begin(), open_b.end(), more_than_b);

    node_f = std::make_shared<Node>(source, std::vector<size_t>(2,0), heuristic(source));
    open_f.push_back(node_f);
    std::push_heap(open_f.begin(), open_f.end(), more_than_f);

    node_b = std::make_shared<Node>(target, std::vector<size_t>(2,0), heuristic_b(target));
    open_b.push_back(node_b);
    std::push_heap(open_b.begin(), open_b.end(), more_than_b);

    while (!open_f.empty() && !open_b.empty())
    {
        // Begin with forward search
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
            std::pop_heap(open_f.begin(), open_f.end(), more_than_f);
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

                open_f.push_back(next);
                std::push_heap(open_f.begin(), open_f.end(), more_than_f);

                closed_f.push_back(node_f);
            }
        }

        // Backward search
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
            std::pop_heap(open_b.begin(), open_b.end(), more_than_b);
            node_b = open_b.back();
            open_b.pop_back();
            num_generation_b +=1;

            // Dominance check
            if ((((1+this->eps[1])*node_b->f_b[1]) >= min_g2_b[source]) ||
                (node_b->g_b[1] >= min_g2_b[node_b->id])) {
                closed_b.push_back(node_b);
                continue;
            }

            min_g2_b[node_b->id] = node_b->g_b[1];
            num_expansion += 1;

            // Find one solution during the search
            if (node_b->id == source) {
                solutions.push_back(node_b);
                log_solution(node_b);
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

                open_f.push_back(next);
                std::push_heap(open_f.begin(), open_f.end(), more_than_f);

                closed_f.push_back(node_f);
            }
        }
    }
}

void BOZStar::update_node(NodePtr node, size_t target, bool is_fw, std::vector<NodePtr>& cur_list,
    const std::vector<NodePtr>& open, SolutionSet &solutions, const Heuristic& heuristic_f) {
    for (const auto& other_node : open) {
        // Initialize a new node representing a solution
        std::vector<size_t> tmp_gf(node->g.size());
        std::vector<size_t> tmp_gb(node->g.size());
        NodePtr parent_f;
        NodePtr parent_b;

        if (is_fw) {
            parent_f = node->parent;
            parent_b = other_node->parent_b;
            for (size_t i = 0; i < tmp_gf.size(); i++) {
                tmp_gf[i] = node->g[i];
                tmp_gb[i] = other_node->g_b[i];
            }
        } else {
            parent_f = other_node->parent;
            parent_b = node->parent_b;
            for (size_t i = 0; i < tmp_gf.size(); i++) {
                tmp_gb[i] = node->g_b[i];
                tmp_gf[i] = other_node->g[i];
            }
        }

        if (node->id == other_node->id) {  // This is a solution
            NodePtr cand_sol = std::make_shared<Node>(target, tmp_gf, tmp_gb, tmp_gb, tmp_gf, 
                parent_f, parent_b);

            bool need_to_add = true;
            std::vector<NodePtr>::iterator it = solutions.begin();
            while (it != solutions.end()) {
                if (is_dominated_dr(cand_sol, *it)) {
                    // candidate node is better than one of the solutions
                    it = solutions.erase(it);
                } else if (is_dominated_dr(*it, cand_sol)) {
                    // one of the solutions is better than candidate node
                    need_to_add = false;
                    break;
                }
            }
            if (need_to_add) {
                solutions.push_back(cand_sol);
            }
        } else {
            // We can update the heuristic of this node
            // Get the front to front heuristic (f2f_h)
            // TODO: A more accurate and fast heuristic coomputation function is needed
            std::vector<size_t> pseudo_h(node->g.size(), 0);
            std::vector<size_t> tmp_f(node->g.size());
            for (size_t i = 0; i < node->g.size(); i++) {
                tmp_f[i] = tmp_gf[i] + tmp_gb[i] + all_pair_path_lbs[node->id][other_node->id][i];
            }
            NodePtr tmp_node = std::make_shared<Node>(target, tmp_f, pseudo_h);

            // Check if this lower bound is dominated by the solution
            bool can_skip = false;
            for (const auto& sol : solutions) {
                if (is_dominated_dr(tmp_node, sol)) {
                    can_skip = true;
                    break;
                }
            }

            if (can_skip) {
                continue;
            } else {
                
            }
        }
    }
}