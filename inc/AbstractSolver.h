#pragma once

#include "Utils/Definitions.h"
#include "Utils/Logger.h"

class AbstractSolver {
protected:
    const AdjacencyMatrix   &adj_matrix;
    // Pair<double>            eps;
    EPS eps;
    double perimeter_dist = 0;
    double perimeter_factor = 0;

    size_t num_sol_f = 0;
    size_t num_sol_b = 0;
    size_t num_expansion = 0;
    size_t num_expansion_f = 0;
    size_t num_expansion_b = 0;
    size_t num_generation= 0;
    size_t num_generation_f = 0;
    size_t num_generation_b = 0;
    size_t num_reinsert = 0;

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
    clock_t runtime_local_check = 0;
    clock_t runtime_global_check = 0;
    size_t branching_factor_f = 0;
    size_t branching_factor_b = 0;
    uint screen;

    virtual std::string get_solver_name() = 0;

    inline size_t get_num_sol_f(){return num_sol_f;}
    inline size_t get_num_sol_b(){return num_sol_b;}
    inline size_t get_num_expansion(){return num_expansion;}
    inline size_t get_num_expansion_f(){return num_expansion_f;}
    inline size_t get_num_expansion_b(){return num_expansion_b;}
    inline size_t get_num_generation(){return num_generation;}
    inline size_t get_num_generation_f(){return num_generation_f;}
    inline size_t get_num_generation_b(){return num_generation_b;}
    inline size_t get_num_reinsert(){return num_reinsert;}
    inline double get_perimeter_dist(void) {return perimeter_dist;}
    inline double get_perimeter_factor(void) {return perimeter_factor;}

    inline void set_perimeter_dist(double val) {perimeter_dist = val;}
    inline void set_perimeter_factor(double val) {perimeter_factor = val;}

    virtual void operator()(size_t source, size_t target, Heuristic &heuristic, 
        SolutionSet &solutions, unsigned int time_limit=UINT_MAX) = 0;

    AbstractSolver(const AdjacencyMatrix &adj_mx, EPS eps, const LoggerPtr logger, uint screen=0): 
        adj_matrix(adj_mx), eps(eps), logger(logger), screen(screen) {}
    virtual ~AbstractSolver(){}
};
