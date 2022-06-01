#ifndef UTILS_DEFINITIONS_H
#define UTILS_DEFINITIONS_H

#include <map>
#include <vector>
#include <array>
#include <list>
#include <iostream>
#include <limits>
#include <functional>
#include <memory>
#include <climits>
#include "boost/heap/pairing_heap.hpp"
#include "boost/heap/priority_queue.hpp"


#ifndef DEBUG
#define DEBUG 1
#endif

using namespace std;

const size_t MAX_COST = std::numeric_limits<size_t>::max();
typedef pair<list<size_t>, vector<size_t>> PathGvalPair;

template<typename T>
using Pair      = std::array<T, 2>;

template<typename T>
std::ostream& operator<<(std::ostream &stream, const Pair<T> pair) {
    stream << "[" << pair[0] << ", " << pair[1] << "]";
    return stream;
}

// template<typename T>
// std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>&b) {
//     assert(a.size() == b.size());
//     std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(a), std::plus<T>());
//     return a;
// }

template<typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>&b) {
    assert(a.size() == b.size());
    std::vector<T> output;
    output.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(output), std::plus<T>());
    return output;
}

template<typename T>
vector<T> get_diff(const vector<T>& a, const vector<T>& b) {
    assert(a.size() == b.size());
    vector<T> output;
    output.reserve(a.size());
    for (size_t t = 0; t < a.size(); t ++) {
        T val = (a[t] < b[t])? b[t] - a[t] : a[t] - b[t];
        output.push_back(val);
    }
    return output;
}

template<typename T>
vector<T> get_comax(const vector<T>& a, const vector<T>& b) {
    assert(a.size() == b.size());
    vector<T> output;
    output.reserve(a.size());
    for (size_t t = 0; t < a.size(); t ++) {
        output.push_back(max(a[t], b[t]));
    }
    return output;
}

using Heuristic = std::function<std::vector<size_t>(size_t)>;

// Structs and classes
struct Edge {
    size_t          source;
    size_t          target;
    std::vector<size_t>    cost;

    Edge(size_t source, size_t target, std::vector<size_t> cost) : source(source), target(target), cost(cost) {}
    Edge inverse() {
        return Edge(this->target, this->source, this->cost);
    }
};
std::ostream& operator<<(std::ostream &stream, const Edge &edge);


// Graph representation as adjacency matrix
class AdjacencyMatrix {
private:
    std::vector<std::vector<Edge>> matrix;
    size_t                         graph_size;
    size_t num_of_objectives = 0;

public:
    AdjacencyMatrix() = default;
    AdjacencyMatrix(size_t graph_size, std::vector<Edge> &edges, bool inverse=false);
    void add(Edge edge);
    size_t size(void) const;
    size_t get_num_of_objectives() const;
    inline size_t get_graph_size(void) const {return graph_size;}
    const std::vector<Edge>& operator[](size_t vertex_id) const;
  
    friend std::ostream& operator<<(std::ostream &stream, const AdjacencyMatrix &adj_matrix);
};


struct Node;
struct PathPair;
struct ApexPathPair;
using NodePtr       = std::shared_ptr<Node>;
using PathPairPtr   = std::shared_ptr<PathPair>;
using ApexPathPairPtr   = std::shared_ptr<ApexPathPair>;
using SolutionSet   = std::vector<NodePtr>;
using PPSolutionSet = std::vector<PathPairPtr>;
using ApexPathSolutionSet = std::vector<ApexPathPairPtr>;
using boost::heap::pairing_heap;
using boost::heap::compare;

using EPS = std::vector<double>;

// TODO: only add a bool variable is_forward to avoid using g_b, h_b, ...
struct Node {
    struct lex_compare {
        bool operator() (const vector<size_t>& a, const vector<size_t>& b) const {
            assert(a->size() == b->size());
            for (size_t i = 0; (i+1) < a.size(); i ++) {
                if (a[i] != b[i]) {
                    return a[i] > b[i];
                }
            }
            return a.back() > b.back();
        }
    };

    size_t id;
    vector<size_t> g;
    vector<size_t> h;
    pairing_heap<vector<size_t>, compare<lex_compare>> other_h;

    vector<size_t> f;
    NodePtr parent;
    list<size_t> path;  // the path that this node represents
    bool is_forward;
    const PathGvalPair* other_path;

    Node (size_t id, vector<size_t> g, vector<size_t> h, 
        NodePtr parent=nullptr, bool is_fwd=true, const PathGvalPair* other_path=nullptr) : 
        id(id), g(g), h(h), f(g+h), parent(parent), is_forward(is_fwd), other_path(other_path) {};
    
    Node (size_t id, vector<size_t> in_g, list<vector<size_t>> in_h=list<vector<size_t>>(), 
        NodePtr parent=nullptr, bool is_fwd=true, const PathGvalPair* other_path=nullptr) : 
        id(id), g(in_g), parent(parent), is_forward(is_fwd), other_path(other_path) {
        if (!in_h.empty()) {
            for (const auto& tmp_h : in_h) {
                other_h.push(tmp_h);
            }
        }
        h = other_h.top();
        f = g + h;
        other_h.pop();
    };

    struct more_than_specific_heurisitic_cost {
        size_t cost_idx;

        more_than_specific_heurisitic_cost(size_t cost_idx) : cost_idx(cost_idx) {};
        bool operator()(const NodePtr &a, const NodePtr &b) const;
    };

    struct more_than_combined_heurisitic {
        double factor;

        more_than_combined_heurisitic(double factor) : factor(factor) {};
        bool operator()(const NodePtr &a, const NodePtr &b) const;
    };


    struct more_than_full_cost {
        bool operator()(const NodePtr &a, const NodePtr &b) const;
    };

    enum LEX_ORDER {LEX0, LEX1};
    struct more_than_lex{
        Node::LEX_ORDER order;
        more_than_lex(Node::LEX_ORDER order) : order(order) {};
        bool operator()(const NodePtr &a, const NodePtr &b) const;
    };


    struct compare_lex1
    {
        bool operator()(const NodePtr n1, const NodePtr n2) const
        {
            if (n1->f[0] != n2->f[0]){
                return n1->f[0] > n2->f[0];
            }
            return n1->f[1] > n2->f[1];
        }
    };
    friend std::ostream& operator<<(std::ostream &stream, const Node &node);

    void set_path(void) {
        if (parent == nullptr) {
            path = list<size_t>({id});
        } else {
            list<size_t> tmp_path = parent->path;
            tmp_path.push_back(id);
            path = tmp_path;
        }
    }

    void set_path(list<size_t> in_path) {
        path = in_path;
    }

    void update_h(void) {
        h = other_h.top();
        other_h.pop();
    }

    void set_hval(list<vector<size_t>> in_h, bool reset=true) {
        if (reset) {
            h.clear();
            other_h.clear();
        } else {
            other_h.push(h);
        }

        for (const auto& tmp_h : in_h) {
            other_h.push(tmp_h);
        }
        update_h();
    }

    void set_hval(vector<size_t> in_h, bool reset=true) {
        if (reset) {
            h.clear();
            other_h.clear();
        } else {
            other_h.push(h);
        }
        other_h.push(in_h);
        update_h();
    }

    void print_h(void) {
        cout << "h val: [";
        for (size_t i = 0; i < h.size(); i++) {
            cout << h[i];
            if (i == h.size()-1) 
                cout << "]" << endl;
            else
                cout << ", ";
        }
        cout << "other_h: " << endl;
        pairing_heap<vector<size_t>, compare<lex_compare>> tmp_h(other_h);
        while(!tmp_h.empty()) {
            vector<size_t> top_h = tmp_h.top();
            tmp_h.pop();
            cout << "[";
            for (size_t i = 0; i < top_h.size(); i++) {
                cout << top_h[i];
                if (i == top_h.size()-1) 
                    cout << "]" << endl;
                else
                    cout << ", ";
            }
        }
        cout << endl;
    }
};


struct PathPair {
    size_t      id;
    NodePtr     top_left;
    NodePtr     bottom_right;
    NodePtr     parent;
    bool        is_active=true;

    PathPair(const NodePtr &top_left, const NodePtr &bottom_right)
        : id(top_left->id), top_left(top_left), bottom_right(bottom_right), parent(top_left->parent) {};

    bool update_nodes_by_merge_if_bounded(const PathPairPtr &other, const Pair<double> eps);
    bool update_nodes_by_merge_if_bounded_keep_track(const PathPairPtr &other, const Pair<double> eps, std::list<NodePtr>& pruned_list);
    bool update_nodes_by_merge_if_bounded2(const PathPairPtr &other, const Pair<double> eps);

    bool if_merge_bounded(const PathPairPtr &other, const Pair<double> eps)  const;


    struct more_than_full_cost {
        bool operator()(const PathPairPtr &a, const PathPairPtr &b) const;
    };

    friend std::ostream& operator<<(std::ostream &stream, const PathPair &pp);
};

enum MergeStrategy {SMALLER_G2, RANDOM, MORE_SLACK, SMALLER_G2_FIRST, REVERSE_LEX};

struct ApexPathPair {
    size_t id; // state of the node
    NodePtr apex;
    NodePtr path_node;
    NodePtr parent;
    Heuristic& h;
    bool is_active=true;

    ApexPathPair(const NodePtr &apex, const NodePtr &path_node, Heuristic& h)
        : id(apex->id), apex(apex), path_node(path_node) , parent(path_node->parent), h(h) {};

    ApexPathPair(const ApexPathPairPtr parent, const Edge& egde);


    bool update_nodes_by_merge_if_bounded(const ApexPathPairPtr &other, const EPS eps, MergeStrategy s=MergeStrategy::SMALLER_G2);
    bool update_apex_by_merge_if_bounded(const NodePtr &other_apex, const EPS eps);

    // bool if_merge_bounded(const ApexPathPairPtr &other, const EP S eps)  const;


    struct more_than_full_cost {
        bool operator()(const ApexPathPairPtr &a, const ApexPathPairPtr &b) const;
    };

    friend std::ostream& operator<<(std::ostream &stream, const ApexPathPair &pp);
};

bool is_bounded(NodePtr apex, NodePtr node,  const EPS eps);
bool is_bounded(NodePtr apex, NodePtr node);
bool is_dominated_dr(NodePtr apex, NodePtr node);
bool is_dominated_dr(NodePtr apex, NodePtr node, const EPS eps);
NodePtr getSource(NodePtr node);


class Interval{
public:
    double eps = 0;
    NodePtr top_left;
    NodePtr bottom_right;
    std::shared_ptr<std::list<NodePtr>> to_expand;

    Interval(){};
    Interval(const NodePtr top_left, const NodePtr bottom_right, std::shared_ptr<std::list<NodePtr>> to_expand);
};

std::ostream& operator<<(std::ostream& os, const Interval& interval);


using IntervalList   = std::vector<Interval>;

typedef boost::heap::priority_queue<NodePtr , boost::heap::compare<Node::compare_lex1> > heap_open_t;

list<PathGvalPair> get_paths(const vector<NodePtr>& in_list);
PathGvalPair combine_path_pair(const PathGvalPair& a, const PathGvalPair& b, const size_t& target);
void floyd_warshell(vector<vector<size_t>>& rst, size_t c_idx, const AdjacencyMatrix& adj_matrix);

#endif //UTILS_DEFINITIONS_H
