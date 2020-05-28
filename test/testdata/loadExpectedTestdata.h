#include <istream>
#include <string>
#include <unordered_map>

#include "MetagraphInterface.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;


auto getExpectedKmers() {
    std::string const testdatapath{TESTDATAPATH};

    // load json
    std::ifstream is(testdatapath + "/expectedKmers.json");
    json j;
    is >> j;

    std::unordered_map<std::string, std::vector<MetagraphInterface::NodeAnnotation>> expectedKmers;
    for (auto&& elem : j.items()) {
        auto kmer = elem.key();
        std::vector<MetagraphInterface::NodeAnnotation> annotations;
        for (auto&& anno : elem.value()) {
            annotations.emplace_back(MetagraphInterface::NodeAnnotation{anno[0].get<std::string>(),
                                                                        anno[1].get<std::string>(),
                                                                        anno[2].get<bool>(),
                                                                        anno[3].get<size_t>()});
        }
        std::sort(annotations.begin(), annotations.end());
        expectedKmers.emplace(kmer, annotations);
    }

    return expectedKmers;
}

struct KmerNeighbours {
    std::vector<std::string> predecessors;
    std::vector<std::string> successors;
};

auto getExpectedNeighbours() {
    std::string const testdatapath{TESTDATAPATH};

    // load json
    std::ifstream is(testdatapath + "/expectedKmerNeighbours.json");
    json j;
    is >> j;

    std::unordered_map<std::string, KmerNeighbours> expectedNeighbours;
    for (auto&& elem : j.items()) {
        auto kmer = elem.key();
        std::vector<std::string> predecessors;
        std::vector<std::string> successors;
        for (auto&& p : elem.value()["prev"]) {
            predecessors.emplace_back(p.get<std::string>());
        }
        for (auto&& s : elem.value()["next"]) {
            successors.emplace_back(s.get<std::string>());
        }
        std::sort(predecessors.begin(), predecessors.end());
        std::sort(successors.begin(), successors.end());
        expectedNeighbours.emplace(kmer, KmerNeighbours{predecessors, successors});
    }

    return expectedNeighbours;
}
