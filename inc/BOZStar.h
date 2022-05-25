#ifndef BI_CRITERIA_BOZ_STAR_H
#define BI_CRITERIA_BOZ_STAR_H

#include <vector>
#include "Utils/Definitions.h"
#include "Utils/Logger.h"
#include "AbstractSolver.h"
#include "BOAStar.h"
#include "ShortestPathHeuristic.h"

class BOZStar: public BOAStar {
private:
    int lookahead_f;
    int lookahead_b;
    std::vector<std::vector<std::vector<size_t>>> all_pair_path_lbs;

protected:
    std::clock_t start_time;
    std::vector<std::pair<std::clock_t, NodePtr>> solution_log;

    size_t num_expansion_f = 0;
    size_t num_generation_f = 0;
    size_t num_expansion_b = 0;
    size_t num_generation_b = 0;

    void init_search(void) {
        num_expansion = 0;
        num_expansion_f = 0;
        num_expansion_b = 0;

        num_generation = 0;
        num_generation_f = 0;
        num_generation_b = 0;
    }

    // is_fw: whether the node is from forward or backward searches
    // cur_list: the list for putting the nodes from the current side
    // open: the open list from the opposite side
    void update_node(NodePtr node, size_t target, bool is_fw, std::vector<NodePtr>& cur_list,
        const std::vector<NodePtr>& open, SolutionSet &solutions, const Heuristic& heuristic_f);

public:
    virtual std::string get_solver_name() {return "BOZ*"; }
    inline int get_look_forward(void) {return lookahead_f;}
    inline int get_look_backward(void) {return lookahead_b;}
    inline void set_look_forward(int val) {lookahead_f = val;}
    inline void set_look_backward(int val) {lookahead_b = val;}

    BOZStar(const AdjacencyMatrix &adj_matrix, Pair<double> eps, const LoggerPtr logger=nullptr, 
        int lh_f=INT_MAX, int lh_b=INT_MAX);

    void operator()(size_t source, size_t target, Heuristic &heuristic, Heuristic &heuristic_b,
        SolutionSet &solutions, AdjacencyMatrix& graph, size_t graph_size, unsigned int time_limit=UINT_MAX);

    std::vector<std::pair<std::clock_t, NodePtr>> get_sol_log(){return solution_log;}
};




#endif //BI_CRITERIA_BOZ_STAR_H
