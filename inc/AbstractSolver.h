#pragma once

#include "Utils/Definitions.h"
#include "Utils/Logger.h"

class AbstractSolver {
protected:
    const AdjacencyMatrix   &adj_matrix;
    // Pair<double>            eps;
    EPS eps;
    size_t perimeter_dist = 0;
    size_t perimeter_factor = 0;

    size_t num_sol_f = 0;
    size_t num_sol_b = 0;
    size_t num_expansion = 0;
    size_t num_expansion_f = 0;
    size_t num_expansion_b = 0;
    size_t num_generation= 0;
    size_t num_generation_f = 0;
    size_t num_generation_b = 0;

    virtual void init_search(){
        num_sol_f = 0;
        num_sol_b = 0;
        num_expansion = 0;
        num_expansion_f = 0;
        num_expansion_b = 0;
        num_generation = 0;
        num_generation_f = 0;
        num_generation_b = 0;
    }

    const LoggerPtr         logger;
    void start_logging(size_t source, size_t target);
    void end_logging(SolutionSet &solutions, bool succ=true);

public:
    clock_t runtime_pre_h = 0;
    clock_t runtime_node_gen = 0;
    clock_t runtime_update_h = 0;
    size_t branching_factor_f = 0;
    size_t branching_factor_b = 0;

    virtual std::string get_solver_name() = 0;

    inline size_t get_num_sol_f(){return num_sol_f;}
    inline size_t get_num_sol_b(){return num_sol_b;}
    inline size_t get_num_expansion(){return num_expansion;}
    inline size_t get_num_expansion_f(){return num_expansion_f;}
    inline size_t get_num_expansion_b(){return num_expansion_b;}
    inline size_t get_num_generation(){return num_generation;}
    inline size_t get_num_generation_f(){return num_generation_f;}
    inline size_t get_num_generation_b(){return num_generation_b;}
    inline size_t get_perimeter_dist(void) {return perimeter_dist;}
    inline size_t get_perimeter_factor(void) {return perimeter_factor;}

    inline void set_perimeter_dist(size_t val) {perimeter_dist = val;}
    inline void set_perimeter_factor(size_t val) {perimeter_factor = val;}

    virtual void operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, unsigned int time_limit=UINT_MAX) = 0;

    AbstractSolver(const AdjacencyMatrix &adj_matrix, EPS eps, const LoggerPtr logger): adj_matrix(adj_matrix), eps(eps), logger(logger) {}
    virtual ~AbstractSolver(){}
};
