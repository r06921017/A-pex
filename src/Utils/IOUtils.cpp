#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include "Utils/IOUtils.h"

void split_string(string in_string, string delimiter, vector<string> &results)
{
    size_t first_delimiter;

    while ((first_delimiter = in_string.find_first_of(delimiter)) != in_string.npos) {
        if (first_delimiter > 0) {
            results.push_back(in_string.substr(0, first_delimiter));
        }
        in_string = in_string.substr(first_delimiter + 1);
    }

    if (in_string.length() > 0) {
        results.push_back(in_string);
    }
}

bool load_gr_files(vector<string> gr_files, vector<Edge> &edges_out, size_t &graph_size){
  size_t max_node_num = 0;
  for (auto gr_file: gr_files){
    ifstream file(gr_file.c_str());
    
    if (file.is_open() == false) {
        cerr << "cannot open the gr file " << gr_file << endl;
        return false;
    }

    string line;
    int idx_edge = 0;
    while (file.eof() == false) {
        getline(file, line);

        if (line == "") {
            break;
        }

        vector<string> decomposed_line;
        split_string(line, " ", decomposed_line);

        string type = decomposed_line[0];
        if ((strcmp(type.c_str(),"c") == 0) || (strcmp(type.c_str(),"p") == 0)) {
            continue; //comment or problem lines, not part of the graph
        }

        if (strcmp(type.c_str(),"a") == 0) { //arc
          if (idx_edge < (int)edges_out.size() - 1){
            if ((stoul(decomposed_line[1]) != edges_out[idx_edge].source) ||
                (stoul(decomposed_line[2]) != edges_out[idx_edge].target)) {
              // arc_sign src dest should be same in both files
              cerr << "file inconsistency" << endl;
              return false;
            }
            edges_out[idx_edge].cost.push_back(stoul(decomposed_line[3]));
          } else {
            Edge e(stoul(decomposed_line[1]),
                   stoul(decomposed_line[2]),
                   {stoul(decomposed_line[3])});
            edges_out.push_back(e);
            max_node_num = max({max_node_num, e.source, e.target});
          }
        }
        idx_edge ++;
    }
    file.close();
  }
  graph_size = max_node_num;
  return true;
}

bool load_gr_files(string gr_file1, string gr_file2, vector<Edge> &edges_out, size_t &graph_size) {
    size_t          max_node_num = 0;
    ifstream   file1(gr_file1.c_str());
    ifstream   file2(gr_file2.c_str());

    if ((file1.is_open() == false) || (file2.is_open() == false)) {
        return false;
    }

    string line1, line2;
    while ((file1.eof() == false) && (file2.eof() == false)) {
        getline(file1, line1);
        getline(file2, line2);

        if ((line1 == "") || (line2 == "")) {
            break;
        }

        vector<string> decomposed_line1, decomposed_line2;
        split_string(line1, " ", decomposed_line1);
        split_string(line2, " ", decomposed_line2);

        string type = decomposed_line1[0];
        if ((strcmp(type.c_str(),"c") == 0) || (strcmp(type.c_str(),"p") == 0)) {
            continue; //comment or problem lines, not part of the graph
        }

        if ((decomposed_line1[0] != decomposed_line2[0]) ||
            (decomposed_line1[1] != decomposed_line2[1]) ||
            (decomposed_line1[2] != decomposed_line2[2])) {
            // arc_sign src dest should be same in both files
            return false;
        }

        if (strcmp(type.c_str(),"a") == 0) { //arc
            Edge e(stoul(decomposed_line1[1]),
                   stoul(decomposed_line1[2]),
                   {stoul(decomposed_line1[3]), stoul(decomposed_line2[3])});
            edges_out.push_back(e);
            max_node_num = max({max_node_num, e.source, e.target});
        }
    }
    graph_size = max_node_num;
    return true;
}

bool load_txt_file(string txt_file, vector<Edge> &edges_out, size_t &graph_size) {
    bool            first_line = true;
    size_t          max_node_num = 0;
    ifstream   file(txt_file.c_str());

    if (file.is_open() == false) {
        return false;
    }

    string line;
    while (file.eof() == false) {
        getline(file, line);

        if (line == "") {
            break;
        }

        vector<string> decomposed_line;
        split_string(line, " ", decomposed_line);

        if (first_line) {
            first_line = false;
            continue;
        }
        Edge e(stoul(decomposed_line[0]),
               stoul(decomposed_line[1]),
               {stoul(decomposed_line[2]), stoul(decomposed_line[3])});
        edges_out.push_back(e);
        max_node_num = max({max_node_num, e.source, e.target});
    }
    graph_size = max_node_num;
    return true;
}


bool load_queries(string query_file, vector<pair<size_t, size_t>> &queries_out) {
    ifstream   file(query_file.c_str());

    if (file.is_open() == false) {
        return false;
    }

    string line;
    while (file.eof() == false) {
        getline(file, line);

        if (line == "") {
            break;
        } else if (line[0] == '#') {
            continue; // Commented out queries
        }

        vector<string> decomposed_line;
        split_string(line, ",", decomposed_line);

        pair<size_t, size_t> query = {stoul(decomposed_line[0]), stoul(decomposed_line[1])};
        queries_out.push_back(query);
    }
    return true;
}
