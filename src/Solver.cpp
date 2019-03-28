#include "Solver.h"

namespace Solver {

std::string solve(const Reader::Graph& graph, const Reader::Parameters& params) {
    std::stringstream results;
    switch(params.solving_mode) {
    case Reader::MODE::NONE:
        // results << "BruteForce:\n";
        // results << brute_force(graph, params);
        // results << "\n\nB&B:\n";
        // results << branch_and_bound(graph, params);
        // results << "\n\nDynamic programming:\n";
        // results << dynamic(graph, params);
        //results << "Tabu:\n";
        //results << tabu(graph, params);
        results << "Genetic:\n";
        results << genetic(graph, params);
        break;
    case Reader::MODE::BRUTE_FORCE:
        results << brute_force(graph, params);
        break;
    case Reader::MODE::BRANCH_AND_BOUND:
        results << branch_and_bound(graph, params);
        break;
    case Reader::MODE::DYNAMIC:
        results << dynamic(graph, params);
        break;
    case Reader::MODE::TABU:
        results << tabu(graph, params);
        break;
    case Reader::MODE::GENETIC:
        results << genetic(graph, params);
        break;
    default:
        throw std::runtime_error("Solving mode not set!");
    }
    return results.str();
}

} // namespace Solver
