#ifndef METAGRAPHINTERFACE_H
#define METAGRAPHINTERFACE_H

#include <functional>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include <tbb/concurrent_queue.h>
#include <graph/annotated_dbg.hpp>
#include <graph/representation/succinct/dbg_succinct.hpp>
#include <annotation/representation/annotation_matrix/static_annotators_def.hpp>
#include <common/utils/file_utils.hpp>

namespace mabl3 {

template <typename V>
void printVector(V const & v) {
    std::cout << "[";
    if (v.size() > 0) {
        std::cout << v.at(0);
        for (size_t i{1}; i < v.size(); ++i) {
            std::cout << ", " << v.at(i);
        }
    }
    std::cout << "]";
}

class MetagraphInterface {
public:
    //! Annotation of a node returned from metagraph
    struct NodeAnnotation {
        std::string sequence;          // fasta header without '> '
        std::vector<size_t> positions; // kmer positions in the sequence
        bool compareSorted(NodeAnnotation const & rhs) const {
            if (sequence == rhs.sequence
                    && positions.size() == rhs.positions.size()) {
                std::vector<size_t> thisPos{};
                std::vector<size_t> rhsPos{};
                thisPos.insert(thisPos.end(), positions.begin(), positions.end());
                rhsPos.insert(rhsPos.end(), rhs.positions.begin(), rhs.positions.end());
                std::sort(thisPos.begin(), thisPos.end());
                std::sort(rhsPos.begin(), rhsPos.end());
                return thisPos == rhsPos;
            }
            return false;
        }
        bool operator==(NodeAnnotation const & rhs) const {
            return sequence == rhs.sequence
                    && positions == rhs.positions;
        }
        bool operator<(NodeAnnotation const & rhs) const {
            return sequence < rhs.sequence
                    || (sequence == rhs.sequence && positions < rhs.positions);
        }
        friend std::ostream& operator<<(std::ostream & out, NodeAnnotation const & annot) {
            out << annot.sequence << ": [";
            if (annot.positions.size()) {
                out << annot.positions.front();
                for (auto it = std::next(annot.positions.begin(), 1); it != annot.positions.end(); ++it) {
                    out << ", " << *it;
                }
            }
            out << "]";
            return out;
        }
    };

    //! Alias for a callback function used in node iteration
    using KmerAnnotationCallback = std::function<void(std::string const &,
                                                      std::vector<NodeAnnotation> const &)>;
    //! Alias for a node identifier
    using NodeID = mtg::graph::DeBruijnGraph::node_index;
    //! Alias for a callback function used in node iteration
    using KmerCallback = std::function<void(std::string const &, NodeID)>;
    //! Alias for a callback function used in node iteration
    using ContigCallback = std::function<void(std::vector<std::string> const &,
                                              std::vector<NodeID> const &)>;
    //! Alias for a callback function used in node iteration
    using ContigAnnotationCallback = std::function<void(std::vector<std::string> const &,
                                                        std::vector<std::vector<NodeAnnotation>> const &)>;
    //! Alias for a callback function used in node iteration
    using ContigAnnotationRawCallback = std::function<void(std::vector<std::string> const &, // kmers
                                                           std::vector<
                                                               std::pair<std::string,        // sequence, [[pos for kmer0], [pos for kmer1], ...]
                                                                         std::vector<SmallVector<uint64_t>>>> const &)>;

    //! c'tor, loads the graph
    /*!
     * \param graph Path to .dbg file
     * \param annotation Path to .row.annodbg file
     */
    MetagraphInterface(std::string const & graphFile,
                       std::string const & annotationFile)
        : graph_{nullptr} {
        // strip file endings from graph paths
        std::string filetype = ".dbg";
        std::string filebase =
            std::equal(filetype.rbegin(), filetype.rend(), graphFile.rbegin())
                ? graphFile.substr(0, graphFile.size() - filetype.size())
                : graphFile;
        std::string annotationFiletype = ".column_coord.annodbg";
        std::string annotationFilebase =
            std::equal(annotationFiletype.rbegin(), annotationFiletype.rend(), annotationFile.rbegin())
                ? annotationFile.substr(0, annotationFile.size() - annotationFiletype.size())
                : annotationFile;
        // load graph and annotation
        auto graph = std::make_shared<mtg::graph::DBGSuccinct>(2); // k doesn't matter, must be at least 2
        if (!graph->load(filebase)) {
            throw std::runtime_error("[ERROR] -- MetagraphInterface -- input file " + graphFile + " corrupted");
        }
        auto annotation = std::make_unique<mtg::annot::ColumnCoordAnnotator>();
        if (!annotation->load(annotationFilebase+".column_coord.annodbg")) {
            throw std::runtime_error("[ERROR] -- MetagraphInterface -- can't load annotations for graph "
                                     + filebase + ".dbg, file "
                                     + annotationFilebase + ".column_coord.annodbg corrupted");
        }
        // store pointer to annotated graph
        graph_ = std::make_unique<mtg::graph::AnnotatedDBG>(graph, std::move(annotation));
        if (!graph_->check_compatibility()) {
            throw std::runtime_error("[ERROR] -- MetagraphInterface -- Graph and Annotation are incompatible.");
        }
        std::cout << numNodes() << " nodes in graph " << graphFile << std::endl;
    }
    //! Get annotation of a node
    std::vector<NodeAnnotation> getAnnotation(NodeID const nodeID) const {
        auto const & coordinates = graph_->get_kmer_coordinates({nodeID}, -1, 0, 0);
        std::vector<NodeAnnotation> coords;
        for (auto&& annotation : coordinates) {
            coords.push_back(parseAnnotation(annotation));
        }
        return coords;
    }
    //! Get annotations of a list of nodes
    std::vector<std::vector<NodeAnnotation>> getAnnotation(std::vector<NodeID> const & nodeIDs) const {
        auto const & coordinates = graph_->get_kmer_coordinates(nodeIDs, -1, 0, 0);
        std::vector<std::vector<NodeAnnotation>> annotations(nodeIDs.size()); // init with same size as nodeIDs
        for (auto&& coord : coordinates) {
            auto const & seq = coord.first;
            for (size_t i{0}; i < coord.second.size(); ++i) {
                if (coord.second[i].size() > 0) {
                    std::vector<size_t> positions{};
                    positions.insert(positions.end(), coord.second[i].begin(), coord.second[i].end());
                    annotations.at(i).emplace_back(NodeAnnotation{seq, positions});
                }
            }
        }
        return annotations;
    }
    //! Get all incoming node IDs
    auto getIncoming(NodeID i) const { return getNeighbors(i, true); }
    //! Get value of k in the graph
    size_t getK() const { return graph_->get_graph().get_k(); }
    //! Given a node index, return the kmer string
    auto getKmer(NodeID i) const { return graph_->get_graph().get_node_sequence(i); }
    //! Get a node from a k-mer string
    NodeID getNode(std::string const & kmer) const { return graph_->get_graph().kmer_to_node(kmer); }
    //! Returns all NodeIDs from a sequence
    std::vector<NodeID> getNodes(std::string const & sequenceName) const {
        std::vector<NodeID> nodes;
        auto callback = [&nodes](NodeID id) { nodes.emplace_back(id); };
        graph_->call_annotated_nodes(sequenceName, callback);
        return nodes;
    }
    //! Returns all NodeIDs that have a sequence and a certain position annotated
    std::vector<NodeID> getNodes(std::string const & sequenceName, size_t pos) const {
        std::vector<NodeID> nodes{};
        for (auto&& nodeID : getNodes(sequenceName)) {
            for (auto&& annot : getAnnotation(nodeID)) {
                if (annot.sequence == sequenceName
                        && std::find(annot.positions.begin(), annot.positions.end(), pos) != annot.positions.end()) {
                    nodes.emplace_back(nodeID);
                }
            }
        }
        return nodes;
    }
    std::vector<NodeID> getNodes(NodeAnnotation const & annot, size_t pos) const {
        return getNodes(annot.sequence, pos);
    }
    //! Get nodes that share an annotation, positions vector needs to be complete!
    auto getNodes(NodeAnnotation const & annot) const {
        std::vector<NodeID> nodes{};
        for (auto&& nodeID : getNodes(annot.sequence)) {
            for (auto&& a : getAnnotation(nodeID)) {
                if (annot.compareSorted(a)) {
                    nodes.emplace_back(nodeID);
                }
            }
        }
        return nodes;
    }
    //! Get all outgoing node IDs
    auto getOutgoing(NodeID i) const { return getNeighbors(i, false); }
    //! Get starting node's IDs
    auto getStartNodes() const {
        std::vector<NodeID> ids;
        graph_->get_graph().call_source_nodes([&ids](NodeID i){ ids.emplace_back(i); });
        return ids;
    }
    //! Direct access to the graph
    auto const & graph() const { return graph_; }
    //! Iterate over all nodes in the graph, calling \c callback for each node. If nthreads > 1, handle any concurrency manually!
    /*! Iterates over all kmers in the graph without filtering */
    void iterateNodes(KmerAnnotationCallback const & callback, size_t nthreads = 1) const {
        iterateGraph([this, &callback](std::string const & kmer, NodeID id){
            callback(kmer, getAnnotation(id));
        }, nthreads);
    }
    //! This overload only returns the kmer string and its NodeID to callback. If nthreads > 1, handle any concurrency manually!
    /*! Iterates over all kmers in the graph without filtering */
    void iterateNodes(KmerCallback const & callback, size_t nthreads = 1) const {
        iterateGraph(callback, nthreads);
    }
    //! This overload returns all kmer strings and the respective NodeIDs of a contig to callback. If nthreads > 1, handle any concurrency manually!
    /*! Iterates over all kmers in the graph without filtering */
    void iterateNodes(ContigAnnotationCallback const & callback, size_t nthreads = 1) const {
        iterateGraph([this, &callback](std::vector<std::string> const & kmers,
                                       std::vector<NodeID> const & ids) {
            callback(kmers, getAnnotation(ids));
        }, nthreads);
    }
    //! This overload returns all kmer strings and the respective raw annotations to callback. If nthreads > 1, handle any concurrency manually!
    /*! Iterates over all kmers in the graph without filtering */
    void iterateNodes(ContigAnnotationRawCallback const & callback, size_t nthreads = 1) const {
        iterateGraph([this, &callback](std::vector<std::string> const & kmers,
                                       std::vector<NodeID> const & ids) {
            callback(kmers, graph_->get_kmer_coordinates(ids, -1, 0, 0));
        }, nthreads);
    }

    //! Get number of nodes in the graph
    size_t numNodes() const { return graph_->get_graph().num_nodes(); }
    //! Get all annotations from metagraph and iterate over them, calling \c callback on each single label
    void parseAllAnnotations(std::function<void(std::string const &)> const & callback) const {
        auto const & annotation = graph_->get_annotation();
        for (auto&& label : annotation.get_all_labels()) { callback(label); }
    }
    //! Push all k-mers on a tbb::concurrent_queue so multiple threads can work with them
    void queueKmers(tbb::concurrent_queue<std::pair<std::string, std::vector<NodeAnnotation>>> & queue) {
        iterateGraph([this, &queue](std::string const & kmer, NodeID id){
            queue.push(std::pair<std::string, std::vector<NodeAnnotation>>(kmer, getAnnotation(id)));
        }, 1); // do not iterate parallel
    }
    auto reconstructSequence(std::string const & sequenceName) const {
        std::unordered_map<size_t, std::string> posToKmer;
        size_t maxPos = 0;
        //size_t k = getK();
        for (auto const & node : getNodes(sequenceName)) {
            auto kmer = getKmer(node);
            for (auto const & annot : getAnnotation(node)) {
                if (annot.sequence == sequenceName) {
                    for (auto const & pos : annot.positions) {
                        posToKmer[pos] = kmer;
                        maxPos = (pos > maxPos) ? pos : maxPos;
                    }
                }
            }
        }
        std::string sequence{};
        size_t i{0};
        while (i <= maxPos) {
            if (posToKmer.find(i) == posToKmer.end() && i == sequence.size()) {
                sequence += "N";
            } else if (posToKmer.find(i) != posToKmer.end()) {
                auto j = sequence.size() - i; // i should never be >sequence.size()
                sequence += posToKmer[i].substr(j);
            }
            ++i;
        }
        return sequence;
    }
    auto reconstructSequence(NodeAnnotation const & annot) const { return reconstructSequence(annot.sequence); }

private:
    //! callback function gets a k-mer and its node ID
    void iterateGraph(KmerCallback const & callback, size_t nthreads) const {
        auto k = getK();
        graph_->get_graph().call_sequences([this,
                                            &callback,
                                            k](std::string const & contig,
                                               std::vector<NodeID> const & ids){
            for (size_t i = 0; i <= contig.size() - k; ++i) {
                callback(contig.substr(i,k), ids[i]);
            }
        },
        nthreads);
    }
    //! callback function gets a vector of k-mers and a vector of the respective node IDs
    void iterateGraph(ContigCallback const & callback, size_t nthreads) const {
        auto k = getK();
        graph_->get_graph().call_sequences([this,
                                            &callback,
                                            k](std::string const & contig,
                                               std::vector<NodeID> const & ids){
            std::vector<std::string> kmers{};
            for (size_t i = 0; i <= contig.size() - k; ++i) {
                kmers.emplace_back(contig.substr(i,k));
            }
            callback(kmers, ids);
        },
        nthreads);
    }
    //! Get all incoming or outcoming \c NodeID s of a certain node
    std::vector<NodeID> getNeighbors(NodeID i, bool incoming) const {
        std::vector<NodeID> neighbours;
        auto callback = [&neighbours](NodeID nodeID) {
            neighbours.emplace_back(nodeID);
        };
        if (incoming) {
            graph_->get_graph().adjacent_incoming_nodes(i, callback);
        } else {
            graph_->get_graph().adjacent_outgoing_nodes(i, callback);
        }
        return neighbours;
    }
    //! Gets a label string from metagraph and returns a NodeAnnotation object
    NodeAnnotation parseAnnotation(std::pair<std::string, std::vector<SmallVector<uint64_t>>> const & annotation) const {
        std::string const & sequence = annotation.first;
        std::vector<size_t> positions{};
        positions.insert(positions.end(), annotation.second.at(0).begin(), annotation.second.at(0).end());
        return NodeAnnotation{sequence, positions};
    }

    //! Pointer to the graph
    std::unique_ptr<mtg::graph::AnnotatedDBG const> graph_;
};

} // namespace mabl3

#endif // METAGRAPHINTERFACE_H
