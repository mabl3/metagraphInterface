#include <chrono>
#include <iostream>
#include <string>

#include "catch2/catch.hpp"
#include "MetagraphInterface.h"
#include <tbb/concurrent_queue.h>
#include "loadExpectedTestdata.h"

TEST_CASE("Metagraph Interface") {
    // load expected data
    auto expectedKmers = getExpectedKmers();
    auto expectedKmerNeighbours = getExpectedNeighbours();
    std::cout << expectedKmers.size() << " kmers in test data" << std::endl;

    // load graph
    std::string const testdatapath{TESTDATAPATH};
    auto graph = MetagraphInterface(testdatapath + "/testdataGraph.dbg",
                                    testdatapath + "/testdataGraph.row.annodbg");
    REQUIRE(graph.getK() == (size_t)39);

    SECTION("Check Annotations") {
        // test parseAllAnnotations
        std::set<MetagraphInterface::NodeAnnotation> observedAnnotationSet;
        std::set<MetagraphInterface::NodeAnnotation> expectedAnnotationSet;
        graph.parseAllAnnotations([&observedAnnotationSet](MetagraphInterface::NodeAnnotation const & annot){
            observedAnnotationSet.emplace(annot);
        });
        // all expeted kmers should be there
        for (auto&& elem : expectedKmers) {
            auto expectedAnnotations = elem.second;
            auto observedAnnotations = graph.getAnnotation(graph.getNode(elem.first));
            std::sort(expectedAnnotations.begin(), expectedAnnotations.end());
            std::sort(observedAnnotations.begin(), observedAnnotations.end());
            REQUIRE(expectedAnnotations == observedAnnotations);
            for (auto&& annot : expectedAnnotations) { expectedAnnotationSet.emplace(annot); }
        }
        REQUIRE(expectedAnnotationSet == observedAnnotationSet);
    }
    SECTION("Check Start Nodes") {
        std::vector<std::string> expectedStart;
        for (auto&& elem : expectedKmerNeighbours) {
            if (elem.second.predecessors.size() == 1
                    && elem.second.predecessors.at(0).at(0) == '$') {
               expectedStart.emplace_back(elem.first);
            }
        }
        std::vector<std::string> observedStart;
        for (auto&& i : graph.getStartNodes()) {
            observedStart.emplace_back(graph.getKmer(i));
        }
        std::sort(expectedStart.begin(), expectedStart.end());
        std::sort(observedStart.begin(), observedStart.end());
        REQUIRE(expectedStart == observedStart);
    }
    SECTION("Check Iteration") {
        // remove softmask and unknown
        auto expectedCleanKmers = expectedKmers;
        for (auto&& elem : expectedKmers) {
            for (char c : elem.first) {
                if (!(c == 'A' || c == 'C' || c == 'G' || c == 'T')) {
                    expectedCleanKmers.erase(elem.first);
                }
            }
        }
        std::unordered_map<std::string, std::vector<MetagraphInterface::NodeAnnotation>> observedKmers;
        auto callback = [&observedKmers](std::string const & kmer,
                                         std::vector<MetagraphInterface::NodeAnnotation> const & occurrences) {
            auto occ = occurrences;
            std::sort(occ.begin(), occ.end());
            observedKmers.emplace(kmer, occ);
        };
        auto ts = std::chrono::system_clock::now();
        graph.iterateNodes(callback, 0, 0, true, true);
        auto tsend = std::chrono::system_clock::now();
        std::cout << "Iteration took " << std::chrono::duration_cast<std::chrono::seconds>(tsend - ts).count() << " s" << std::endl;
        REQUIRE(expectedCleanKmers == observedKmers);
    }
    SECTION("Check Fast Iteration") {
        std::unordered_map<std::string, std::vector<MetagraphInterface::NodeAnnotation>> observedKmers;
        auto callback = [&observedKmers](std::string const & kmer,
                                         std::vector<MetagraphInterface::NodeAnnotation> const & occurrences) {
            auto occ = occurrences;
            std::sort(occ.begin(), occ.end());
            observedKmers.emplace(kmer, occ);
        };
        auto ts = std::chrono::system_clock::now();
        graph.iterateNodes(callback, 0, 0, false, false);
        auto tsend = std::chrono::system_clock::now();
        std::cout << "Fast iteration took " << std::chrono::duration_cast<std::chrono::seconds>(tsend - ts).count() << " s" << std::endl;
        REQUIRE(expectedKmers == observedKmers);
    }
    SECTION("Check Queueing") {
        tbb::concurrent_queue<std::pair<std::string, std::vector<MetagraphInterface::NodeAnnotation>>> queue;
        graph.queueKmers(queue);
        std::unordered_map<std::string, std::vector<MetagraphInterface::NodeAnnotation>> observedKmers;
        std::pair<std::string, std::vector<MetagraphInterface::NodeAnnotation>> elem;
        while (queue.try_pop(elem)) {
            auto occ = elem.second;
            std::sort(occ.begin(), occ.end());
            observedKmers.emplace(elem.first, occ);
        }
        REQUIRE(expectedKmers == observedKmers);
    }
    SECTION("Check Neighbours") {
        for (auto&& elem : expectedKmerNeighbours) {
            std::vector<std::string> observedIncoming;
            std::vector<std::string> observedOutgoing;
            for (auto&& i : graph.getIncoming(graph.getNode(elem.first))) {
                observedIncoming.emplace_back(graph.getKmer(i));
            }
            for (auto&& i : graph.getOutgoing(graph.getNode(elem.first))) {
                observedOutgoing.emplace_back(graph.getKmer(i));
            }
            std::sort(observedIncoming.begin(), observedIncoming.end());
            std::sort(observedOutgoing.begin(), observedOutgoing.end());
            REQUIRE(elem.second.predecessors == observedIncoming);
            REQUIRE(elem.second.successors == observedOutgoing);
        }
    }
}
