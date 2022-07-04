#include <iostream>
#include <memory>
#include <time.h>
#include <fstream>

#include "ShortestPathHeuristic.h"
#include "Utils/Definitions.h"
#include "Utils/IOUtils.h"
#include "Utils/Logger.h"
#include "BOAStar.h"
#include "PPA.h"
#include "SingleCriteria.h"
#include "ApexSearch.h"
#include "NAMOA.h"
#include "BOPS.h"

#include <boost/program_options.hpp>
#include<boost/tokenizer.hpp>

const string resource_path = "resources/";
const string output_path = "output/";
const MergeStrategy DEFAULT_MERGE_STRATEGY = MergeStrategy::SMALLER_G2;
string alg_variant = "";


// Simple example to demonstarte the usage of the algorithm
void run_algo(size_t graph_size, AdjacencyMatrix& graph, AdjacencyMatrix&inv_graph, 
    size_t source, size_t target, ofstream& output, string algorithm, MergeStrategy ms, 
    LoggerPtr logger, double eps, double peri_factor, uint b_mode, uint time_limit, uint screen) {

    // Compute heuristic
    if (screen > DEBUG_LOG_EXPANSION) {
        cout << "Source: " << source << ", Target: " << target << endl;
        cout << "Start computing heuristic... ";
    }
    clock_t pre_start = clock();
    ShortestPathHeuristic sp_heuristic(target, graph_size, inv_graph);
    double runtime_pre_process = (double)(clock() - pre_start) / CLOCKS_PER_SEC;
    if (screen > DEBUG_LOG_EXPANSION) cout << "done!" << endl;

    using placeholders::_1;
    Heuristic heuristic = bind( &ShortestPathHeuristic::operator(), sp_heuristic, _1);

    SolutionSet solutions;
    auto runtime = clock();
    auto start =clock();

    unique_ptr<AbstractSolver> solver;
    if (algorithm == "PPA"){
        Pair<double> eps_pair({eps, eps});
        solver = make_unique<PPA>(graph, eps_pair, logger);
    } else if (algorithm == "BOA"){
        Pair<double> eps_pair({eps, eps});
        solver = make_unique<BOAStar>(graph, eps_pair, logger, screen);
    } else if (algorithm == "NAMOAdr"){
        EPS eps_vec (graph.get_num_of_objectives(), eps);
        solver = make_unique<NAMOAdr>(graph, eps_vec, logger);
        // ((ApexSearch*)solver.get())->set_merge_strategy(ms);
    } else if (algorithm == "Apex"){
        EPS eps_vec (graph.get_num_of_objectives(), eps);
        solver = make_unique<ApexSearch>(graph, eps_vec, logger);
        ((ApexSearch*)solver.get())->set_merge_strategy(ms);
    } else if (algorithm == "BOPS") {
        pre_start = clock();
        ShortestPathHeuristic sp_heuristic_b(source, graph_size, inv_graph);
        Heuristic heuristic_b = bind( &ShortestPathHeuristic::operator(), sp_heuristic_b, _1);
        runtime_pre_process += (double)(clock() - pre_start) / CLOCKS_PER_SEC;
        Pair<double> eps_pair({eps, eps});
        solver = make_unique<BOPS>(graph, eps_pair, heuristic, heuristic_b, 
            peri_factor, b_mode, logger, screen);
    } else {
        cerr << "unknown solver name" << endl;
        exit(-1);
    }

    (*solver)(source, target, heuristic, solutions, time_limit);
    runtime = clock() - start;

    string sol_name = solver->get_solver_name();
    size_t num_gen = solver->get_num_generation();
    size_t num_gen_f = solver->get_num_generation_f();
    size_t num_gen_b = solver->get_num_generation_b();
    size_t num_exp = solver->get_num_expansion();
    size_t num_exp_f = solver->get_num_expansion_f();
    size_t num_exp_b = solver->get_num_expansion_b();
    size_t num_reinsert = solver->get_num_reinsert();
    double runtime_pre_h = (double) solver->runtime_pre_h / CLOCKS_PER_SEC;
    double runtime_node_gen = (double) solver->runtime_node_gen / CLOCKS_PER_SEC;
    double runtime_update_h = (double) solver->runtime_update_h / CLOCKS_PER_SEC;
    double runtime_global_check = (double) solver->runtime_global_check / CLOCKS_PER_SEC;
    double runtime_local_check = (double) solver->runtime_local_check / CLOCKS_PER_SEC;
    size_t bf_f = solver->branching_factor_f;  // forward branching factor
    size_t bf_b = solver->branching_factor_b;  // backward branching factor

    if (runtime / CLOCKS_PER_SEC > time_limit) cout << "Timeout :( ";
    else cout << "Success :) ";
    cout << "Runtime: " <<  (double) runtime / CLOCKS_PER_SEC <<
        ", #exp: " << num_exp << ", #gen: " << num_gen << ", #sol: " << solutions.size() << endl;
    if (screen > DEBUG_LOG_EXPANSION)
        for (auto sol: solutions) cout << *sol << endl;

    output << (double) runtime / CLOCKS_PER_SEC << "," << source << "," << target << "," << 
        sol_name << "," << eps << "," << peri_factor << "," << solutions.size() << "," <<
        num_gen << "," << num_gen_f << "," << num_gen_b << "," <<
        num_exp << "," << num_exp_f << "," << num_exp_b << "," <<
        runtime_pre_process << "," << runtime_pre_h << "," << runtime_node_gen << "," << 
        runtime_update_h << "," << (runtime_global_check + runtime_local_check) << "," << 
        runtime_global_check << "," << runtime_local_check << "," <<
        bf_f << "," << bf_b << "," << num_reinsert;

    output << endl;
    if (screen > DEBUG_LOG_EXPANSION) cout << "-----End Single Example-----" << endl;
}

void single_run_map(size_t graph_size, vector<Edge> & edges, size_t source, size_t target, 
    string output_file, string algorithm, MergeStrategy ms, LoggerPtr logger, double eps, 
    double peri_factor, uint b_mode, uint time_limit, uint screen) {

    AdjacencyMatrix graph(graph_size, edges);
    AdjacencyMatrix inv_graph(graph_size, edges, true);

    ifstream infile(output_path + output_file);
    bool exist = infile.good();
    infile.close();
    if (!exist) {
        ofstream addHeads(output_path + output_file);
        addHeads << "runtime,source,target," << 
            "algorithm,eps,factor,#solutions," <<
            "#generation,#forward generation,#backward generation," << 
            "#expansion,#forward expansion,#backward expansion," << 
            "runtime preprocess h,runtime build perimeter,runtime node generation," << 
            "runtime update h,runtime dominance check," << 
            "runtime global dominance check,runtime local dominance check," <<
            "#forward branching factor,#backward branching factor,#reinsert";
        addHeads << endl;
    }

    ofstream stats;
    stats.open(output_path + output_file, fstream::app);

    run_algo(graph_size, graph, inv_graph, source, target, stats, algorithm, ms, logger, eps, 
        peri_factor, b_mode, time_limit, screen);
 }

void run_query(size_t graph_size, vector<Edge> & edges, string query_file, string output_file, 
    string algorithm, MergeStrategy ms, LoggerPtr logger, double eps, 
    double peri_factor, uint b_mode, uint time_limit, uint screen) {

    ifstream infile(output_path + output_file);
    bool exist = infile.good();
    infile.close();
    if (!exist) {
        ofstream addHeads(output_path + output_file);
        addHeads << "runtime,source,target," << 
            "algorithm,eps,factor,#solutions," <<
            "#generation,#forward generation,#backward generation," << 
            "#expansion,#forward expansion,#backward expansion," << 
            "runtime preprocess h,runtime build perimeter,runtime node generation," << 
            "runtime update h,runtime dominance check," << 
            "runtime global dominance check,runtime local dominance check," <<
            "#forward branching factor,#backward branching factor,#reinsert";
        addHeads << endl;
    }

    ofstream stats;
    stats.open(output_path + output_file, fstream::app);

    vector<pair<size_t, size_t>> queries;
    if (load_queries(query_file, queries) == false) {
        cout << "Failed to load queries file" << endl;
        return;
    }

    // Build graphs
    AdjacencyMatrix graph(graph_size, edges);
    AdjacencyMatrix inv_graph(graph_size, edges, true);

    size_t query_count = 0;
    for (auto iter = queries.begin(); iter != queries.end(); ++iter) {
        query_count++;
        cout << "Started Query: " << query_count << "/" << queries.size() << endl;
        size_t source = iter->first;
        size_t target = iter->second;

        run_algo(graph_size, graph, inv_graph, source, target, stats, algorithm, ms, logger, 
            eps, peri_factor, b_mode, time_limit, screen);
    }

}

int main(int argc, char** argv){
    namespace po = boost::program_options;

    vector<string> objective_files;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("start,s", po::value<int>()->default_value(-1), "start location")
        ("goal,g", po::value<int>()->default_value(-1), "goal location")
        ("query,q", po::value<string>()->default_value(""), "number of agents")
        ("map,m",po::value<vector<string> >(&objective_files)->multitoken(), "files for edge weight")
        ("eps,e", po::value<double>()->default_value(0), "approximation factor")
        ("factor,f", po::value<double>()->default_value(4.0), "perimeter factor")
        ("backward", po::value<uint>()->default_value(0), "backward solver mode (0:BOD, 1:BOA*)")
        ("merge", po::value<string>()->default_value(""), "strategy for merging apex node pair: SMALLER_G2, RANDOM or MORE_SLACK")
        ("algorithm,a", po::value<string>()->default_value("Apex"), "solvers (BOA, PPA or Apex search)")
        ("cutoffTime,t", po::value<uint>()->default_value(300), "cutoff time (seconds)")
        ("output,o", po::value<string>()->required(), "Name of the output file")
        ("logging_file", po::value<string>()->default_value(""), "logging file" )
        ("debug", po::value<uint>()->default_value(0), "debug log shown on the terminal" )
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        cout << desc << endl;
        return 1;
    }

    po::notify(vm);
    srand((int)time(0));

    uint screen = vm["debug"].as<uint>();

    if (vm["query"].as<string>() != ""){
        cout << "Query file: " << vm["query"].as<string>() << endl;
        if (vm["start"].as<int>() != -1 || vm["goal"].as<int>() != -1){
            cerr << "query file and start/goal cannot be given at the same time !" << endl;
            return -1;
        }
    }
    
    LoggerPtr logger = nullptr;

    if (vm["logging_file"].as<string>() != ""){
        logger = new Logger(vm["logging_file"].as<string>());
    }

    // Load files
    size_t graph_size;
    vector<Edge> edges;

    if (screen > DEBUG_LOG_EXPANSION)
        for (auto file:objective_files)
            cout << file << endl;

    if (load_gr_files(objective_files, edges, graph_size) == false) {
        cout << "Failed to load gr files" << endl;
        return -1;
    }

    if (screen > DEBUG_LOG_EXPANSION) cout << "Graph Size: " << graph_size << endl;

    // Build graphs
    MergeStrategy ms = DEFAULT_MERGE_STRATEGY;
    alg_variant = vm["merge"].as<string>();

    if (vm["merge"].as<string>() != "" && vm["algorithm"].as<string>()!= "Apex"){
        alg_variant = "";
        cout << "WARNING: merge strategy with non-apex search" << endl;
    }else if(vm["merge"].as<string>() == "SMALLER_G2"){
        ms = MergeStrategy::SMALLER_G2;
    }else if(vm["merge"].as<string>() == "SMALLER_G2_FIRST"){
        ms = MergeStrategy::SMALLER_G2_FIRST;
    }else if(vm["merge"].as<string>() == "RANDOM"){
        ms = MergeStrategy::RANDOM;
    }else if(vm["merge"].as<string>() == "MORE_SLACK"){
        ms = MergeStrategy::MORE_SLACK;
    }else if(vm["merge"].as<string>() == "REVERSE_LEX"){
        ms = MergeStrategy::REVERSE_LEX;
    }else if (screen > DEBUG_LOG_EXPANSION){
        cout << "unknown merge strategy" << endl;
    }

    if (vm["query"].as<string>() != ""){
        run_query(graph_size, edges, vm["query"].as<string>(), 
            vm["output"].as<string>(), vm["algorithm"].as<string>(), ms, logger, 
            vm["eps"].as<double>(), vm["factor"].as<double>(), vm["backward"].as<uint>(), 
            vm["cutoffTime"].as<uint>(), screen);
    } else {
        single_run_map(graph_size, edges, vm["start"].as<int>(), vm["goal"].as<int>(), 
            vm["output"].as<string>(), vm["algorithm"].as<string>(), ms, logger, 
            vm["eps"].as<double>(), vm["factor"].as<double>(), vm["backward"].as<uint>(), 
            vm["cutoffTime"].as<uint>(), screen);
    }

    delete(logger);

    return 0;
}
