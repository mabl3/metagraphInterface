#include <chrono>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "catch2/catch.hpp"
#include "MetagraphInterface.h"
#include <tbb/concurrent_queue.h>
#include "loadExpectedTestdata.h"

/* Output specializations for easier debugging in case of errors */
template <typename V>
std::ostream& operator<<(std::ostream & out, std::vector<V> vector) {
    out << "[";
    if (vector.size()) {
        out << vector.front();
        for (auto it = std::next(vector.begin(), 1); it != vector.end(); ++it) {
            out << ", " << *it;
        }
    }
    out << "]";
    return out;
}

template <typename K, typename V>
std::ostream& operator<<(std::ostream & out, std::unordered_map<K, V> map) {
    out << "{";
    if (map.size()) {
        out << map.begin()->first << ": " << map.begin()->second;
        for (auto it = std::next(map.begin(), 1); it != map.end(); ++it) {
            out << ",\n    " << it->first << ": " << it->second;
        }
    }
    out << "}";
    return out;
}

template <typename T>
inline bool mapEqual(T const & h1, T const & h2) {
    auto ret = true;
    for (auto& elem : h1) {
        if (h2.find(elem.first) == h2.end()) {
            std::cout << "'" << elem.first << "' not found in h2" << std::endl;
            ret = false;
        } else if (h2.at(elem.first) != elem.second) {
            std::cout << "For '" << elem.first << "', " << h2.at(elem.first) << " (h2) != " << elem.second << " (h1)" << std::endl;
            ret = false;
        }
    }
    for (auto& elem : h2) {
        if (h1.find(elem.first) == h1.end()) {
            std::cout << "'" << elem.first << "' not found in h1" << std::endl;
            ret = false;
        } else if (h1.at(elem.first) != elem.second) {
            std::cout << "For '" << elem.first << "', " << h1.at(elem.first) << " (h1) != " << elem.second << " (h2)" << std::endl;
            ret = false;
        }
    }
    return ret;
}

/* Custom Hash Functions for NodeAnnotation */
inline void combineHash(size_t & seed, size_t hash) {
    seed ^= hash + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

struct AnnotationHash {
    size_t operator()(mabl3::MetagraphInterface::NodeAnnotation const & a) const {
        size_t seed = 0;
        combineHash(seed, std::hash<std::string>{}(a.sequence));
        for (auto&& pos : a.positions) {
            combineHash(seed, std::hash<size_t>{}(pos));
        }
        return seed;
    }
};

//! Implements operator() to check if two Link s are equal, ignores span and score
struct AnnotationEqual {
    bool operator()(mabl3::MetagraphInterface::NodeAnnotation const & lhs,
                    mabl3::MetagraphInterface::NodeAnnotation const & rhs) const {
        return lhs == rhs;
    }
};



TEST_CASE("Metagraph Interface") {
    std::cout << "[Test Case: Metagraph Interface]" << std::endl;
    // load graph
    std::string const testdatapath{TESTDATAPATH};
    auto graph = mabl3::MetagraphInterface(testdatapath + "/testdataGraph.dbg",
                                           testdatapath + "/testdataGraph.column_coord.annodbg");

    // load expected data
    auto expectedKmers = getExpectedKmers();
    auto expectedKmerNeighbours = getExpectedNeighbours();
    std::cout << expectedKmers.size() << " kmers in test data" << std::endl;
    REQUIRE(graph.getK() == (size_t)12);
    SECTION("Check Iteration") {
        std::cout << "[Section: Check Iteration]" << std::endl;
        std::unordered_map<std::string, std::vector<mabl3::MetagraphInterface::NodeAnnotation>> observedKmers;
        auto callback = [&observedKmers](std::string const & kmer,
                                         std::vector<mabl3::MetagraphInterface::NodeAnnotation> const & occurrences) {
            auto occ = occurrences;
            std::sort(occ.begin(), occ.end());
            observedKmers.emplace(kmer, occ);
        };
        auto ts = std::chrono::system_clock::now();
        graph.iterateNodes(callback);
        auto tsend = std::chrono::system_clock::now();
        std::cout << "Iteration took " << std::chrono::duration_cast<std::chrono::seconds>(tsend - ts).count() << " s" << std::endl;
        REQUIRE(mapEqual(expectedKmers, observedKmers));
        // --- Only kmer iteration ---
        std::unordered_set<std::string> expectedKmerSet;
        for (auto&& elem : expectedKmers) { expectedKmerSet.emplace(elem.first); }
        std::unordered_set<std::string> observedKmerSet;
        auto kmerCallback = [&observedKmerSet](std::string const & kmer, mabl3::MetagraphInterface::NodeID nodeID) {
            (void) nodeID;
            observedKmerSet.emplace(kmer);
        };
        auto ts2 = std::chrono::system_clock::now();
        graph.iterateNodes(kmerCallback);
        auto ts2end = std::chrono::system_clock::now();
        std::cout << "Iteration (only kmers) took " << std::chrono::duration_cast<std::chrono::seconds>(ts2end - ts2).count() << " s" << std::endl;
        REQUIRE(expectedKmerSet == observedKmerSet);
    }
    SECTION("Check Queueing") {
        std::cout << "[Section: Check Queueing]" << std::endl;
        tbb::concurrent_queue<std::pair<std::string, std::vector<mabl3::MetagraphInterface::NodeAnnotation>>> queue;
        graph.queueKmers(queue);
        std::unordered_map<std::string, std::vector<mabl3::MetagraphInterface::NodeAnnotation>> observedKmers;
        std::pair<std::string, std::vector<mabl3::MetagraphInterface::NodeAnnotation>> elem;
        while (queue.try_pop(elem)) {
            auto occ = elem.second;
            std::sort(occ.begin(), occ.end());
            observedKmers.emplace(elem.first, occ);
        }
        REQUIRE(expectedKmers == observedKmers);
    }
    SECTION("Check Neighbours") {
        std::cout << "[Section: Check Neighbours]" << std::endl;
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
    SECTION("Check Start Nodes") {
        std::cout << "[Section: Check Start Nodes]" << std::endl;
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
    SECTION("Check Nodes from Annotation") {
        // *** TODO: test getNode(s) functions ***
        std::cout << "[Section: Check Nodes from Annotation]" << std::endl;
        std::unordered_map<std::string,
                           std::unordered_set<mabl3::MetagraphInterface::NodeID>> seqToNodes{};
        std::unordered_map<mabl3::MetagraphInterface::NodeAnnotation,
                           std::unordered_set<mabl3::MetagraphInterface::NodeID>,
                           AnnotationHash, AnnotationEqual> annoToNodes{};
        std::unordered_map<mabl3::MetagraphInterface::NodeAnnotation,
                           mabl3::MetagraphInterface::NodeID,
                           AnnotationHash, AnnotationEqual> posToNode{};

        // Fill maps with expected query results
        auto callbackFill = [&graph,
                             &seqToNodes,
                             &annoToNodes,
                             &posToNode](std::string const & kmer,
                                         mabl3::MetagraphInterface::NodeID id) {
            REQUIRE(graph.getNode(kmer) == id);
            auto annotations = graph.getAnnotation(id);
            for (auto&& annotation : annotations) {
                seqToNodes[annotation.sequence].insert(id);
                annoToNodes[annotation].insert(id);
                for (auto&& pos : annotation.positions) {
                    mabl3::MetagraphInterface::NodeAnnotation posAnno{annotation.sequence,
                                                                      std::vector<mabl3::MetagraphInterface::NodeID>{pos}};
                    REQUIRE(posToNode.find(posAnno) == posToNode.end()); // each position should be uniquely annotated
                    posToNode[posAnno] = id;
                }
            }
        };
        graph.iterateNodes(callbackFill);

        // check getNodes(seq) with sequenceToNodes
        for (auto&& elem : seqToNodes) {
            auto ids = graph.getNodes(elem.first);
            std::unordered_set<mabl3::MetagraphInterface::NodeID> idSet{};
            idSet.insert(ids.begin(), ids.end());
            REQUIRE(idSet == elem.second);
        }

        // iterating over graph nodes, getting kmer and annotations per node, testing whether reverse querying works
        size_t i = 0;
        auto callback = [&graph, &i,
                         &annoToNodes,
                         &posToNode](std::string const & kmer,
                                       std::vector<mabl3::MetagraphInterface::NodeAnnotation> const & annotations) {
            if (i < 1000) { // checking all nodes takes way too long
                auto trueID = graph.getNode(kmer);
                for (auto&& annotation : annotations) {
                    for (auto&& pos : annotation.positions) {
                        // test if seq+pos lead to correct id
                        auto ids = graph.getNodes(annotation, pos);
                        REQUIRE(ids.size() == 1);
                        REQUIRE(ids.at(0) == trueID);
                    }
                    // get ids that share the exact annotation
                    auto ids = graph.getNodes(annotation);
                    std::unordered_set<mabl3::MetagraphInterface::NodeID> obsIDs{};
                    obsIDs.insert(ids.begin(), ids.end());
                    REQUIRE(obsIDs == annoToNodes.at(annotation));
                    REQUIRE(obsIDs.find(trueID) != obsIDs.end());
                    /*for (auto&& id : ids) {
                        std::unordered_set<std::string> seqs;
                        for (auto&& anno : graph.getAnnotation(id)) { seqs.insert(anno.sequence); }
                        REQUIRE(std::find(seqs.begin(), seqs.end(), annotation.sequence) != seqs.end());  // make sure all ids share the current sequence
                    }*/
                }
                ++i;
            }
        };
        graph.iterateNodes(callback);
    }
    SECTION("Check Annotations") {
        std::cout << "[Section: Check Annotations]" << std::endl;
        std::set<std::string> observedLabelSet;
        std::set<std::string> expectedLabelSet;
        for (auto&& elem : expectedKmers) {
            for (auto&& annot : elem.second) {
                expectedLabelSet.insert(annot.sequence);
            }
        }
        graph.parseAllAnnotations([&observedLabelSet](std::string const & label){
            observedLabelSet.emplace(label);
        });
        REQUIRE(expectedLabelSet == observedLabelSet);
    }
}
