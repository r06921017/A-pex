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

using namespace std;

const string resource_path = "resources/";
const string output_path = "output/";
const MergeStrategy DEFAULT_MERGE_STRATEGY = MergeStrategy::SMALLER_G2;
string alg_variant = "";


// Simple example to demonstarte the usage of the algorithm
void single_run_map(size_t graph_size, AdjacencyMatrix& graph, AdjacencyMatrix&inv_graph, 
    size_t source, size_t target, ofstream& output, string algorithm, MergeStrategy ms, 
    LoggerPtr logger, double eps, size_t peri_factor, unsigned int time_limit) {

    // Compute heuristic
    cout << "Start Computing Heuristic" << endl;
    clock_t pre_start = clock();
    ShortestPathHeuristic sp_heuristic(target, graph_size, inv_graph);
    double runtime_pre_process = (double)(clock() - pre_start) / CLOCKS_PER_SEC;
    cout << "Finish Computing Heuristic\n" << endl;

    using std::placeholders::_1;
    Heuristic heuristic = bind( &ShortestPathHeuristic::operator(), sp_heuristic, _1);

    SolutionSet solutions;
    auto runtime = std::clock();
    auto start =std::clock();

    std::unique_ptr<AbstractSolver> solver;
    if (algorithm == "PPA"){
        Pair<double> eps_pair({eps, eps});
        solver = make_unique<PPA>(graph, eps_pair, logger);
    } else if (algorithm == "BOA"){
        Pair<double> eps_pair({eps, eps});
        solver = make_unique<BOAStar>(graph, eps_pair, logger);
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
        cout << "heuristic_b(target).front(): " << heuristic_b(target).front() << endl;
        runtime_pre_process += (double)(clock() - pre_start) / CLOCKS_PER_SEC;
        Pair<double> eps_pair({eps, eps});
        solver = make_unique<BOPS>(graph, eps_pair, heuristic, heuristic_b, logger, peri_factor);
    } else {
        cerr << "unknown solver name" << endl;
        exit(-1);
    }

    (*solver)(source, target, heuristic, solutions, time_limit);
    runtime = clock() - start;

    cout << "Node expansion: " << solver->get_num_expansion() << endl;
    cout << "Runtime: " <<  (double) runtime / CLOCKS_PER_SEC << endl;
    cout << "#solutions: " << solutions.size() << endl;
    for (auto sol: solutions) cout << *sol << endl;

    string sol_name = solver->get_solver_name();
    size_t num_gen = solver->get_num_generation();
    size_t num_gen_f = solver->get_num_generation_f();
    size_t num_gen_b = solver->get_num_generation_b();
    size_t num_exp = solver->get_num_expansion();
    size_t num_exp_f = solver->get_num_expansion_f();
    size_t num_exp_b = solver->get_num_expansion_b();
    double runtime_pre_h = (double) solver->runtime_pre_h / CLOCKS_PER_SEC;
    double runtime_node_gen = (double) solver->runtime_node_gen / CLOCKS_PER_SEC;
    double runtime_update_h = (double) solver->runtime_update_h / CLOCKS_PER_SEC;
    size_t bf_f = solver->branching_factor_f;  // forward branching factor
    size_t bf_b = solver->branching_factor_b;  // backward branching factor

    output << (double) runtime / CLOCKS_PER_SEC << "," << source << "," << target << "," << 
        sol_name << "," << eps << "," << peri_factor << "," << solutions.size() << "," <<
        num_gen << "," << num_gen_f << "," << num_gen_b << "," <<
        num_exp << "," << num_exp_f << "," << num_exp_b << "," <<
        runtime_pre_h << "," << runtime_node_gen << "," << runtime_update_h << "," <<
        bf_f << "," << bf_b;

    output << endl;
    cout << "-----End Single Example-----" << endl;
}

void single_run_map(size_t graph_size, vector<Edge> & edges, size_t source, size_t target, 
    string output_file, string algorithm, MergeStrategy ms, LoggerPtr logger, double eps, 
    size_t peri_factor, int time_limit) {

    AdjacencyMatrix graph(graph_size, edges);
    AdjacencyMatrix inv_graph(graph_size, edges, true);

    ifstream infile(output_path + output_file);
    bool exist = infile.good();
    infile.close();
    if (!exist) {
        ofstream addHeads(output_path + output_file);
        addHeads << "runtime,source,target,algorithm,eps,factor,#solutions," <<
            "#generation,#forward generation,#backward generation," << 
            "#expansion,#forward expansion,#backward expansion," << 
            "runtime preprocess h,runtime build perimeter,runtime node generation,runtime update h," << 
            "#forward branching factor,#backward branching factor";
        addHeads << endl;
    }

    ofstream stats;
    stats.open(output_path + output_file, fstream::app);

    single_run_map(graph_size, graph, inv_graph, source, target, stats, algorithm, ms, logger, eps, 
        peri_factor, time_limit);
 }

void run_query(size_t graph_size, vector<Edge> & edges, string query_file, string output_file, 
    string algorithm, MergeStrategy ms, LoggerPtr logger, double eps, size_t peri_factor, int time_limit) {
    ofstream stats;
    stats.open(output_path + output_file, std::fstream::app);

    std::vector<std::pair<size_t, size_t>> queries;
    if (load_queries(query_file, queries) == false) {
        std::cout << "Failed to load queries file" << std::endl;
        return;
    }

    // Build graphs
    AdjacencyMatrix graph(graph_size, edges);
    AdjacencyMatrix inv_graph(graph_size, edges, true);

    size_t query_count = 0;
    for (auto iter = queries.begin(); iter != queries.end(); ++iter) {

        query_count++;
        std::cout << "Started Query: " << query_count << "/" << queries.size() << std::endl;
        size_t source = iter->first;
        size_t target = iter->second;

        single_run_map(graph_size, graph, inv_graph, source, target, stats, algorithm, ms, logger, 
            eps, peri_factor, time_limit);
    }

}

int main(int argc, char** argv){
    namespace po = boost::program_options;

    std::vector<string> objective_files;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("start,s", po::value<int>()->default_value(-1), "start location")
        ("goal,g", po::value<int>()->default_value(-1), "goal location")
        ("query,q", po::value<std::string>()->default_value(""), "number of agents")
        ("map,m",po::value< std::vector<string> >(&objective_files)->multitoken(), "files for edge weight")
        ("eps,e", po::value<double>()->default_value(0), "approximation factor")
        ("factor,f", po::value<size_t>()->default_value(4), "perimeter factor")
        ("merge", po::value<std::string>()->default_value(""), "strategy for merging apex node pair: SMALLER_G2, RANDOM or MORE_SLACK")
        ("algorithm,a", po::value<std::string>()->default_value("Apex"), "solvers (BOA, PPA or Apex search)")
        ("cutoffTime,t", po::value<uint>()->default_value(300), "cutoff time (seconds)")
        ("output,o", po::value<std::string>()->required(), "Name of the output file")
        ("logging_file", po::value<std::string>()->default_value(""), "logging file" )
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    po::notify(vm);
    srand((int)time(0));

    if (vm["query"].as<std::string>() != ""){
        if (vm["start"].as<int>() != -1 || vm["goal"].as<int>() != -1){
            std::cerr << "query file and start/goal cannot be given at the same time !" << std::endl;
            return -1;
        }
    }
    
    LoggerPtr logger = nullptr;

    if (vm["logging_file"].as<std::string>() != ""){
        logger = new Logger(vm["logging_file"].as<std::string>());
    }

    // Load files
    size_t graph_size;
    std::vector<Edge> edges;

    for (auto file:objective_files){
        std::cout << file << std::endl;
    }


    if (load_gr_files(objective_files, edges, graph_size) == false) {
        std::cout << "Failed to load gr files" << std::endl;
        return -1;
    }

    std::cout << "Graph Size: " << graph_size << std::endl;

    // Build graphs
    MergeStrategy ms = DEFAULT_MERGE_STRATEGY;
    alg_variant = vm["merge"].as<std::string>();

    if (vm["merge"].as<std::string>() != "" && vm["algorithm"].as<std::string>()!= "Apex"){
        alg_variant = "";
        std::cout << "WARNING: merge strategy with non-apex search" << std::endl;
    }else if(vm["merge"].as<std::string>() == "SMALLER_G2"){
        ms = MergeStrategy::SMALLER_G2;
    }else if(vm["merge"].as<std::string>() == "SMALLER_G2_FIRST"){
        ms = MergeStrategy::SMALLER_G2_FIRST;
    }else if(vm["merge"].as<std::string>() == "RANDOM"){
        ms = MergeStrategy::RANDOM;
    }else if(vm["merge"].as<std::string>() == "MORE_SLACK"){
        ms = MergeStrategy::MORE_SLACK;
    }else if(vm["merge"].as<std::string>() == "REVERSE_LEX"){
        ms = MergeStrategy::REVERSE_LEX;
    }else{
        std::cerr << "unknown merge strategy" << std::endl;
    }


    if (vm["query"].as<std::string>() != ""){
        run_query(graph_size, edges, vm["query"].as<std::string>(), vm["output"].as<std::string>(), 
            vm["algorithm"].as<std::string>(), ms, logger, vm["eps"].as<double>(), 
            vm["factor"].as<uint>(), vm["cutoffTime"].as<uint>());
    } else{
        cout << "factor: " << vm["factor"].as<size_t>() << endl;
        single_run_map(graph_size, edges, vm["start"].as<int>(), vm["goal"].as<int>(), 
            vm["output"].as<std::string>(), vm["algorithm"].as<std::string>(), ms, logger, 
            vm["eps"].as<double>(), vm["factor"].as<size_t>(), vm["cutoffTime"].as<uint>());
    }

    delete(logger);

    return 0;
}
