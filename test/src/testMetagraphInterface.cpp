#include <iostream>
#include <tuple>
#include <unordered_map>

#include "catch2/catch.hpp"
#include "MetagraphInterface.h"
#include "loadExpectedTestdata.h"

TEST_CASE("Metagraph Interface") {
    auto expectedKmers = getExpectedKmers();
    auto expectedKmerNeighbours = getExpectedNeighbours();

    for (auto&& elem : expectedKmerNeighbours) {
        std::cout << elem.first << "\n  predecessors:" << std::endl;
        for (auto&& p : elem.second.predecessors) {
            std::cout << "\t" << p << std::endl;
        }
        std::cout << "  successors:" << std::endl;
        for (auto&& s : elem.second.successors) {
            std::cout << "\t" << s << std::endl;
        }
    }
}
