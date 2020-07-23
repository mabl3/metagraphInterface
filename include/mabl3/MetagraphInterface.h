#ifndef METAGRAPHINTERFACE_H
#define METAGRAPHINTERFACE_H

#include <functional>
#include <iostream>
#include <unordered_set>

#include <tbb/concurrent_queue.h>
#include <graph/annotated_dbg.hpp>
#include <graph/representation/succinct/dbg_succinct.hpp>
#include <annotation/representation/row_compressed/annotate_row_compressed.hpp>
#include <common/utils/file_utils.hpp>

class MetagraphInterface {
public:
    //! Annotation of a node returned from metagraph
    struct NodeAnnotation {
        std::string genome;     // filename (e.g. hg38.fasta)
        std::string sequence;   // fasta header without '> '
        bool reverse_strand;
        size_t bin_idx;         // kmer position in the sequence, possibly rounded
        bool operator==(NodeAnnotation const & rhs) const {
            return genome == rhs.genome
                    && sequence == rhs.sequence
                    && reverse_strand == rhs.reverse_strand
                    && bin_idx == rhs.bin_idx;
        }
        bool operator<(NodeAnnotation const & rhs) const {
            return genome < rhs.genome
                    || (genome == rhs.genome && sequence < rhs.sequence)
                    || (genome == rhs.genome && sequence == rhs.sequence && bin_idx < rhs.bin_idx)
                    || (genome == rhs.genome && sequence == rhs.sequence && bin_idx == rhs.bin_idx
                        && reverse_strand < rhs.reverse_strand);
        }
    };
    struct RedundantNodeAnnotation : NodeAnnotation {
        RedundantNodeAnnotation(std::string const & _genome,
                                std::string const & _sequence,
                                bool _reverse_strand,
                                size_t _bin_idx,
                                std::string const & _annotationString)
            : NodeAnnotation{_genome, _sequence, _reverse_strand, _bin_idx},
              annotationString{_annotationString} {}
        RedundantNodeAnnotation(NodeAnnotation const & _annotation,
                                std::string const & _annotationString)
            : NodeAnnotation{_annotation}, annotationString{_annotationString} {}
        std::string annotationString;
    };

    //! Alias for a callback function used in node iteration
    using KMerCallback = std::function<void(std::string const &,
                                            std::vector<NodeAnnotation> const &)>;
    //! Alias for a node identifier
    using NodeID = DeBruijnGraph::node_index;

    //! c'tor, loads the graph
    /*!
     * \param graph Path to .dbg file
     * \param annotation Path to .row.annodbg file
     */
    MetagraphInterface(std::string const & graphFile,
                       std::string const & annotationFile)
        : graph_{nullptr} {
        // strip file endings from paths
        std::string filetype = ".dbg";
        std::string filebase =
            std::equal(filetype.rbegin(), filetype.rend(), graphFile.rbegin())
                ? graphFile.substr(0, graphFile.size() - filetype.size())
                : graphFile;
        std::string annotationFiletype = ".row.annodbg";
        std::string annotationFilebase =
            std::equal(annotationFiletype.rbegin(), annotationFiletype.rend(), annotationFile.rbegin())
                ? annotationFile.substr(0, annotationFile.size() - annotationFiletype.size())
                : annotationFile;
        // load graph and annotation
        auto graph = std::make_shared<DBGSuccinct>(2); // k doesn't matter, must be at least 2
        if (!graph->load(filebase)) {
            throw std::runtime_error("[ERROR] -- MetagraphInterface -- input file " + graphFile + " corrupted");
        }
        auto annotation = std::make_unique<annotate::RowCompressed<std::string>>(0, false);
        if (!annotation->load(annotationFilebase)) {
            throw std::runtime_error("[ERROR] -- MetagraphInterface -- can't load annotations for graph "
                                     + filebase + ".dbg, file "
                                     + annotationFilebase + ".row.annodbg corrupted");
        }
        // store pointer to annotated graph
        graph_ = std::make_unique<AnnotatedDBG>(graph,
                                                std::unique_ptr<annotate::MultiLabelEncoded<std::string>>(annotation.release()));
        if (!graph_->check_compatibility()) {
            throw std::runtime_error("[ERROR] -- MetagraphInterface -- Graph and Annotation are incompatible.");
        }
        std::cout << numNodes() << " nodes in graph " << graphFile << std::endl;
    }
    //! Get annotation of a node
    std::vector<NodeAnnotation> getAnnotation(NodeID const nodeID) const {
        auto const & labels = graph_->get_labels(nodeID);
        std::vector<NodeAnnotation> bins;
        for (std::string const & label : labels) {
            bins.push_back(parseLabel(label));
        }
        return bins;
    }
    //! Get redundant annotation of a node
    std::vector<RedundantNodeAnnotation> getRedundantAnnotation(NodeID const nodeID) const {
        auto const & labels = graph_->get_labels(nodeID);
        std::vector<RedundantNodeAnnotation> bins;
        for (std::string const & label : labels) {
            bins.push_back(RedundantNodeAnnotation{parseLabel(label), label});
        }
        return bins;
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
        return getNodes(mergeLabel(annot));
    }
    //! Get nodes that share an annotation
    auto getNodes(RedundantNodeAnnotation const & annot) const {
        return getNodes(annot.annotationString);
    }
    //! Returns all NodeIDs that share an annotation
    std::vector<NodeID> getNodes(std::string const & annotationstring) const {
        std::vector<NodeID> nodes;
        auto callback = [&nodes](NodeID id) { nodes.emplace_back(id); };
        graph_->call_annotated_nodes(annotationstring, callback);
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
    void iterateNodes(KMerCallback const & callback,
                          size_t minNumGenomes = 0,
                          size_t maxKmerCount = 0,
                          bool onlyACGT = false,
                          bool onlyUC = false) const {
        if (minNumGenomes || maxKmerCount || onlyACGT || onlyUC) {
            std::unordered_set<std::string> genomes;    // keep track for minNumGenomes
            iterateGraph([this,
                              &callback,
                              &genomes,
                              &minNumGenomes,
                              &maxKmerCount,
                              &onlyACGT,
                              &onlyUC](std::string const & kmer, NodeID id){
                genomes.clear();
                if (onlyACGT && kmer.find('N') != std::string::npos) { return; }
                if (onlyUC) {
                    for (char c : kmer) {
                        // a lower case letter detected
                        if (c > 90) { return; }
                    }
                }
                auto bins = getAnnotation(id);
                for (auto const & bin : bins) {
                    genomes.insert(bin.genome);
                }
                // if all requirements fulfilled, execute the callback function on this node
                if (genomes.size() >= minNumGenomes
                        && (maxKmerCount == 0 || bins.size() - genomes.size() < maxKmerCount)) {
                    callback(kmer, bins);
                }
            });
        } else {
            // no need for checking, just return each kmer to this methods callback
            iterateGraph([this,
                              &callback](std::string const & kmer, NodeID id){
                callback(kmer, getAnnotation(id));
            });
        }
    }
    //! Get number of nodes in the graph
    size_t numNodes() const { return graph_->get_graph().num_nodes(); }
    //! Get all annotations from metagraph and iterate over them, calling \c callback on each single label
    void parseAllAnnotations(std::function<void(NodeAnnotation const & annot)> const & callback) const {
        auto const & annotation = graph_->get_annotation();
        for (auto&& label : annotation.get_all_labels()) {
            callback(parseLabel(label));
        }
    }
    //! Push all k-mers on a tbb::concurrent_queue so multiple threads can work with them
    void queueKmers(tbb::concurrent_queue<std::pair<std::string, std::vector<NodeAnnotation>>> & queue) const {
        iterateGraph([this, &queue](std::string const & kmer, NodeID id){
            queue.push(std::pair<std::string, std::vector<NodeAnnotation>>(kmer, getAnnotation(id)));
        });
    }

    friend std::ostream& operator<<(std::ostream & out, NodeAnnotation const & bin) {
        out << bin.genome << ", " << bin.sequence << ", " << bin.bin_idx << ", " << bin.reverse_strand;
        return out;
    }

private:
    //! callback function gets a k-mer and its node ID
    void iterateGraph(std::function<void(std::string const &, NodeID)> const & callback) const {
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
            graph_->get_graph().adjacent_incoming_nodes(i,callback);
        } else {
            graph_->get_graph().adjacent_outgoing_nodes(i, callback);
        }
        return neighbours;
    }
    //! Gets a NodeAnnotation and creates a metagraph label string
    std::string mergeLabel(NodeAnnotation const & annotation) const {
        return annotation.genome + "\1" + annotation.sequence + "\1"
                + std::to_string(annotation.reverse_strand) + "\1"
                + std::to_string(annotation.bin_idx);
    }
    //! Gets a label string from metagraph and returns a NodeAnnotation object
    NodeAnnotation parseLabel(std::string const & label) const {
        const auto & parts = utils::split_string(label, "\1");
        if (parts.size() != 4) {
            throw std::runtime_error("[ERROR] -- MetagraphInterface::parseLabel -- invalid label " + label);
        }
        std::string const & genome = parts[0];
        std::string const & sequence = parts[1];
        bool reverse = std::stoi(parts[2]) != 0;
        size_t bin_idx = static_cast<size_t>(std::stoi(parts[3]));
        return NodeAnnotation{genome, sequence, reverse, bin_idx};
    }

    //! Pointer to the graph
    std::unique_ptr<AnnotatedDBG> graph_;
};

#endif // METAGRAPHINTERFACE_H
