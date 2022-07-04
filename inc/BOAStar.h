#ifndef BI_CRITERIA_BOA_STAR_H
#define BI_CRITERIA_BOA_STAR_H

#include <vector>
#include "Utils/Definitions.h"
#include "Utils/Logger.h"
#include "AbstractSolver.h"

class BOAStar: public AbstractSolver {
protected:
    clock_t start_time;

    vector<pair<clock_t, NodePtr>> solution_log;
    void log_solution(NodePtr);

public:
    virtual string get_solver_name() {return "BOA*"; }

    BOAStar(const AdjacencyMatrix &adj_matrix, Pair<double> eps, 
        const LoggerPtr logger=nullptr, uint screen=0);

    void operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, 
        uint time_limit=UINT_MAX);

    vector<pair<clock_t, NodePtr>> get_sol_log(){return solution_log;}
};




#endif //BI_CRITERIA_BOA_STAR_H
