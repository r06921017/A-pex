#ifndef BI_CRITERIA_BOZ_STAR_H
#define BI_CRITERIA_BOZ_STAR_H

#include <vector>
#include "Utils/Definitions.h"
#include "Utils/Logger.h"
#include "AbstractSolver.h"
#include "BOAStar.h"
#include "ShortestPathHeuristic.h"

class BOPS: public BOAStar {
private:
    size_t lookahead_f;
    size_t lookahead_b;
    std::vector<std::vector<std::vector<size_t>>> all_pair_lbs;

protected:
    std::clock_t start_time;
    std::vector<std::pair<std::clock_t, NodePtr>> solution_log;

    size_t num_expansion_f = 0;
    size_t num_generation_f = 0;
    size_t num_expansion_b = 0;
    size_t num_generation_b = 0;
    size_t num_better_sol = 0;
    size_t num_redundent_sol = 0;

    void init_search(void) {
        num_expansion = 0;
        num_expansion_f = 0;
        num_expansion_b = 0;

        num_generation = 0;
        num_generation_f = 0;
        num_generation_b = 0;

        num_better_sol = 0;
        num_redundent_sol = 0;
    }

    // is_fw: whether the node is from forward or backward searches
    // cur_list: the list for putting the nodes from the current side
    // open: the open list from the opposite side
    void update_open(vector<NodePtr>& open, const vector<NodePtr>& other_open, 
        const vector<NodePtr>& closed, const vector<NodePtr>& other_closed,
        size_t target, SolutionSet &solutions);

public:
    virtual std::string get_solver_name() {return "BOPS"; }
    inline size_t get_look_forward(void) {return lookahead_f;}
    inline size_t get_look_backward(void) {return lookahead_b;}
    inline void set_look_forward(size_t val) {lookahead_f = val;}
    inline void set_look_backward(size_t val) {lookahead_b = val;}

    BOPS(const AdjacencyMatrix &adj_matrix, Pair<double> eps, const LoggerPtr logger=nullptr, 
        size_t lh_f=SIZE_MAX, size_t lh_b=SIZE_MAX);

    void operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet & solutions, 
        unsigned int time_limit);

    // void operator()(size_t source, size_t target, SolutionSet &solutions, AdjacencyMatrix& graph, 
    //     size_t graph_size, unsigned int time_limit=UINT_MAX);

    std::vector<std::pair<std::clock_t, NodePtr>> get_sol_log(){return solution_log;}
};




#endif //BI_CRITERIA_BOZ_STAR_H
