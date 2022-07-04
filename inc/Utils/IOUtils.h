#ifndef UTILS_IO_UTILS_H
#define UTILS_IO_UTILS_H

#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include "Utils/Definitions.h"

using namespace std;

bool load_gr_files(string gr_file1, string gr_file2, vector<Edge> &edges, size_t &graph_size);
bool load_gr_files(vector<string> gr_files, vector<Edge> &edges, size_t &graph_size);
bool load_queries(string query_file, vector<pair<size_t, size_t>> &queries_out);

#endif //UTILS_IO_UTILS_H
