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
            // check redundant annotations
            auto observedRedundantAnnotation = graph.getRedundantAnnotation(graph.getNode(elem.first));
            for (size_t i = 0; i < observedAnnotations.size(); ++i) {
                auto const & annostr = observedRedundantAnnotation.at(i).annotationString;
                auto parts = utils::split_string(annostr, "\1");
                MetagraphInterface::NodeAnnotation reconstructed{parts.at(0),
                                                                 parts.at(1),
                                                                 static_cast<bool>(std::stoi(parts.at(2))),
                                                                 static_cast<size_t>(std::stoi(parts.at(3)))};
                REQUIRE(observedAnnotations.at(i).bin_idx == observedRedundantAnnotation.at(i).bin_idx);
                REQUIRE(observedAnnotations.at(i).genome == observedRedundantAnnotation.at(i).genome);
                REQUIRE(observedAnnotations.at(i).reverse_strand == observedRedundantAnnotation.at(i).reverse_strand);
                REQUIRE(observedAnnotations.at(i).sequence == observedRedundantAnnotation.at(i).sequence);
                REQUIRE(observedAnnotations.at(i) == reconstructed);
            }
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
        // --- Only kmer iteration ---
        std::unordered_set<std::string> expectedKmerSet;
        for (auto&& elem : expectedKmers) { expectedKmerSet.emplace(elem.first); }
        std::unordered_set<std::string> observedKmerSet;
        auto kmerCallback = [&observedKmerSet](std::string const & kmer, MetagraphInterface::NodeID nodeID) {
            (void) nodeID;
            observedKmerSet.emplace(kmer);
        };
        auto ts2 = std::chrono::system_clock::now();
        graph.iterateNodes(kmerCallback);
        auto ts2end = std::chrono::system_clock::now();
        std::cout << "Fast iteration (only kmers) took " << std::chrono::duration_cast<std::chrono::seconds>(ts2end - ts2).count() << " s" << std::endl;
        REQUIRE(expectedKmerSet == observedKmerSet);
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
    SECTION("Check Nodes from Annotation") {
        size_t i = 0;
        auto callback = [&graph, &i](std::string const & kmer,
                                     std::vector<MetagraphInterface::NodeAnnotation> const & annotations) {
            if (i < 10000) { // checking all nodes takes way too long
                auto trueID = graph.getNode(kmer);
                for (auto&& annotation : annotations) {
                    auto ids = graph.getNodes(annotation);
                    REQUIRE(std::find(ids.begin(), ids.end(), trueID) != ids.end());    // find trueID in ids
                    for (auto&& id : ids) {
                        auto annos = graph.getAnnotation(id);
                        REQUIRE(std::find(annos.begin(), annos.end(), annotation) != annos.end());  // make sure all ids share the current annotation
                    }
                }
                ++i;
            }
        };
        graph.iterateNodes(callback);
    }
}
