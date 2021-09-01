#ifndef METAGRAPHINTERFACE_H
#define METAGRAPHINTERFACE_H

#include <functional>
#include <iostream>
#include <unordered_set>

#include <tbb/concurrent_queue.h>
#include <graph/annotated_dbg.hpp>
#include <graph/representation/succinct/dbg_succinct.hpp>
#include <annotation/representation/annotation_matrix/static_annotators_def.hpp>
#include <common/utils/file_utils.hpp>

namespace mabl3 {

class MetagraphInterface {
public:
    //! Annotation of a node returned from metagraph
    struct NodeAnnotation {
        std::string sequence;          // fasta header without '> '
        std::vector<size_t> positions; // kmer positions in the sequence
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
    //! Get all incoming node IDs
    auto getIncoming(NodeID i) const { return getNeighbors(i, true); }
    //! Get value of k in the graph
    size_t getK() const { return graph_->get_graph().get_k(); }
    //! Given a node index, return the kmer string
    auto getKmer(NodeID i) const { return graph_->get_graph().get_node_sequence(i); }
    //! Get a node from a k-mer string
    NodeID getNode(std::string const & kmer) const { return graph_->get_graph().kmer_to_node(kmer); }
    //! Get nodes that share an annotation
    auto getNodes(NodeAnnotation const & annot) const {
        //return getNodes(mergeLabel(annot));
        return getNodes(annot.sequence);
    }
    //! Returns all NodeIDs from a sequence
    std::vector<NodeID> getNodes(std::string const & sequence) const {
        std::vector<NodeID> nodes;
        auto callback = [&nodes](NodeID id) { nodes.emplace_back(id); };
        graph_->call_annotated_nodes(sequence, callback);
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
    //! Iterate over all nodes in the graph, calling \c callback for each node
    /*! Iterates over all kmers in the graph without filtering */
    void iterateNodes(KmerAnnotationCallback const & callback) const {
        iterateGraph([this, &callback](std::string const & kmer, NodeID id){
            callback(kmer, getAnnotation(id));
        });
    }
    //! This overload only returns the kmer string and its NodeID to callback
    /*! Iterates over all kmers in the graph without filtering */
    void iterateNodes(KmerCallback const & callback) const {
        iterateGraph(callback);
    }
    //! Get number of nodes in the graph
    size_t numNodes() const { return graph_->get_graph().num_nodes(); }
    //! Get all annotations from metagraph and iterate over them, calling \c callback on each single label
    void parseAllAnnotations(std::function<void(std::string const & sequence)> const & callback) const {
        auto const & annotation = graph_->get_annotation();
        for (auto&& label : annotation.get_all_labels()) { callback(label); }
    }
    //! Push all k-mers on a tbb::concurrent_queue so multiple threads can work with them
    void queueKmers(tbb::concurrent_queue<std::pair<std::string, std::vector<NodeAnnotation>>> & queue) const {
        iterateGraph([this, &queue](std::string const & kmer, NodeID id){
            queue.push(std::pair<std::string, std::vector<NodeAnnotation>>(kmer, getAnnotation(id)));
        });
    }

private:
    //! callback function gets a k-mer and its node ID
    void iterateGraph(KmerCallback const & callback) const {
        auto k = getK();
        graph_->get_graph().call_sequences([this,
                                            &callback,
                                            k](std::string const & contig,
                                               std::vector<NodeID> const & ids){
            for (size_t i = 0; i <= contig.size() - k; ++i) {
                callback(contig.substr(i,k), ids[i]);
            }
        });
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
