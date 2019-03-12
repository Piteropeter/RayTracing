#pragma once
#include "Reader.h"
#include "BruteForce.h"
#include "BranchAndBound.h"
#include "Dynamic.h"
#include "Tabu.h"
#include "Genetic.h"
#include <string>

namespace Solver {

std::string solve(const Reader::Graph& graph, const Reader::Parameters& params);

}