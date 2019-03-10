#include "Reader.h"

namespace Reader {

void mirror_matrix(std::vector<std::vector<double>>& matrix) {
    for(auto i = 0u; i < matrix.size(); i++) {
        for(auto j = 0u; j < matrix.size(); j++) {
            matrix[i][j] = matrix[j][i];
        }
    }
}

void insert_element(const std::string& element, std::istream& file, Graph& graph) {
    std::string buffer;
    if(element == "N:") {
        file >> buffer;
        graph.nodes_count = std::stoul(buffer);
    } else if(element == "M:") {
        file >> buffer;
        graph.modules_count = std::stoul(buffer);
    } else if(element == "Symmetric:") {
        file >> buffer;
        if(buffer == "true")
            graph.symmetric = true;
        else if(buffer == "false")
            graph.symmetric = false;
        else
            throw std::runtime_error("\nReader error: Symmetric not set!\n");
    } else if(element == "Triangle:") {
        file >> buffer;
        if(buffer == "true")
            graph.triangle = true;
        else if(buffer == "false")
            graph.triangle = false;
        else
            throw std::runtime_error("\nReader error: Triangle not set!\n");
    } else {
        graph.modules.resize(graph.modules_count);

        auto size = std::stoul(element);
        graph.modules[0].resize(size);
        for(auto i = 0u; i < size; i++) {
            file >> buffer;
            graph.modules[0][i] = std::stoul(buffer);
        }

        for(auto i = 1u; i < graph.modules_count; i++) {
            file >> buffer;
            auto size = std::stoul(buffer);
            graph.modules[i].resize(size);
            for(auto j = 0u; j < size; j++) {
                file >> buffer;
                graph.modules[i][j] = std::stoul(buffer);
            }
        }

        graph.distance_matrix.resize(graph.nodes_count);
        for(auto i = 0u; i < graph.nodes_count; i++) {
            graph.distance_matrix[i].resize(graph.nodes_count);
            for(auto j = 0u; j < graph.nodes_count; j++) {
                file >> buffer;
                if(i == j) {
                    graph.distance_matrix[i][j] = MAX_VALUE;
                    continue;
                }
                graph.distance_matrix[i][j] = std::stoul(buffer);
            }
        }
    }
}

std::string extract_name(const std::string& path) {
    std::size_t begin = path.rfind('\\');
    std::size_t end = path.find(".txt");

    return path.substr(begin + 1, end - begin - 1);
}

Graph read_graph(const std::string& path) {
    std::ifstream file;
    std::string buffer;
    Graph graph;

    try {
        file.open(path);
        file >> buffer;

        while(!file.eof()) {
            insert_element(buffer, file, graph);
            file >> buffer;
        }

        file.close();

        graph.name = extract_name(path);
    } catch(std::exception& e) {
        std::cout << std::endl << e.what();
    }

    return graph;
}
} // namespace Reader
