#ifndef BI_CRITERIA_BOPS_H
#define BI_CRITERIA_BOPS_H

#include <vector>
#include "Utils/Definitions.h"
#include "Utils/Logger.h"
#include "AbstractSolver.h"
#include "BOAStar.h"
#include "ShortestPathHeuristic.h"

class BOPS: public BOAStar {
private:
    Heuristic heuristic_f;
    Heuristic heuristic_b;
    vector<NodePtr> perimeter;

protected:
    std::clock_t start_time;
    std::vector<std::pair<std::clock_t, NodePtr>> solution_log;

    void init_search(void) {
        num_expansion = 0;
        num_expansion_f = 0;
        num_expansion_b = 0;

        num_generation = 0;
        num_generation_f = 0;
        num_generation_b = 0;
    }

    vector<size_t> get_diff_heuristic(size_t loc1, size_t loc2);
    void reinsert(NodePtr node, vector<NodePtr>& open, vector<NodePtr>& closed, 
        const Node::more_than_full_cost& more_than);
    void print_gcl(const vector<vector_heap>& in_gcl, int id);

public:
    virtual string get_solver_name() {return "BOPS"; }

    BOPS(const AdjacencyMatrix &adj_matrix, Pair<double> eps, Heuristic &h_f, Heuristic &h_b, 
        const LoggerPtr logger=nullptr, size_t perimeter_factor=2);

    void operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet & solutions, 
        unsigned int time_limit);

    vector<pair<clock_t, NodePtr>> get_sol_log(){return solution_log;}
};

#endif //BI_CRITERIA_BOPS_H
