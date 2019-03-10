#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

namespace Reader {

constexpr std::size_t MAX_VALUE = 0xffffffff;

enum class MODE { NONE, BRUTE_FORCE, BRANCH_AND_BOUND, DYNAMIC, TABU, GENETIC };

struct Parameters {
    bool is_path_set = false;
    bool benchmark_mode = false;
    bool help_mode = false;
    bool file_mode = false;
    bool diversification = false;
    double time = 1.0;
    double crossing_chance = 0.9;
    double mutation_chance = 0.05;
    MODE solving_mode = MODE::NONE;
    std::string path = "";
    std::size_t seed = 2018;
    std::size_t population_size = 1000;
    std::size_t tabu_list_size = 100;
    std::size_t benchmark_repetitions = 1;
    std::size_t algorithm_iterations = 100;
};

struct Graph {
    std::string name = "";
    std::size_t nodes_count = 0;
    std::size_t modules_count = 0;
    bool symmetric = false;
    bool triangle = false;
    std::vector<std::vector<std::size_t>> modules;
    std::vector<std::vector<std::size_t>> distance_matrix;
};

void insert_element(const std::string& element, std::istream& file, Graph& graph);
std::string extract_name(const std::string& path);
Graph read_graph(const std::string& name);

} // namespace Reader
