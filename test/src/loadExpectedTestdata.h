#include <filesystem>
#include <istream>
#include <string>
#include <unordered_map>

#include "MetagraphInterface.h"
#include "nlohmann/json.hpp"

namespace fs = std::filesystem;
using json = nlohmann::json;



auto getExpectedKmers() {
    std::string const testdatapath{TESTDATAPATH};

    // load json
    std::ifstream is(testdatapath + "/expectedKmers.json");
    json j;
    is >> j;

    std::unordered_map<std::string, std::vector<mabl3::MetagraphInterface::NodeAnnotation>> expectedKmers;
    for (auto&& elem : j.items()) {
        auto kmer = elem.key();
        std::vector<mabl3::MetagraphInterface::NodeAnnotation> annotations;
        for (auto&& anno : elem.value()) {
            annotations.emplace_back(mabl3::MetagraphInterface::NodeAnnotation{anno[0].get<std::string>(),
                                                                               anno[1].get<std::vector<size_t>>()});
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



auto getSequences() {
    auto chomp = [](std::string& line) {
        if (line.size() > 0) {
            // chomp \n and \r if present ( https://en.wikipedia.org/wiki/Newline#Representation )
            while (line.back() == '\n' || line.back() == '\r') {
                line.pop_back();
            }
        }
    };

    std::unordered_map<std::string, std::string> headToSeq;
    fs::path const testdatapath{TESTDATAPATH};
    for (auto const & entry : fs::directory_iterator{testdatapath}) {
        auto const & f = entry.path();
        if (f.extension() == ".fa" || f.extension() == ".fasta") {
            // Try to open file (non binary mode)
            std::ifstream inputStream(f);
            if (!inputStream.is_open()) {
                throw std::runtime_error("[ERROR] -- getSequences -- Failed to open " + f.string());
            }
            // If file is opened successfully, try to read contents
            std::string line;
            std::string currentHead;
            size_t lineCount{0};
            // getline returns inputStream and as https://en.cppreference.com/w/cpp/io/basic_ios/operator_bool,
            //     body only gets executed if getline was successful (i.e. no fail or read after eof)
            while ( std::getline(inputStream, line) ) {
                ++lineCount;
                chomp(line);    // strip newline sequences from line
                if (line.size() > 0) {  // skip empty lines
                    if (line.at(0) == '>') { // new sequence
                        line.erase(line.begin());   // remove ">"
                        currentHead = line;
                        // new (empty) sequence entry
                        headToSeq.insert({currentHead, ""});
                    } else if (line.at(0) == ';') { // skip comment line
                        continue;
                    } else {    // sequence line
                        if (currentHead.empty()) {
                            auto msg = "[ERROR] -- getSequences -- Encountered non-header, non-comment line"
                                       " that does not belong to any sequence header in '" + f.string() +
                                       "' in line " + std::to_string(lineCount);
                            throw std::runtime_error(msg);
                        }
                        headToSeq.at(currentHead).append(line);    // create contiguous sequence from respective entries
                    }
                }
            }
            if ( !inputStream.eof() && inputStream.fail() ) {
                auto msg = "[ERROR] -- getSequences -- Failed to read " + f.string();
                throw std::runtime_error(msg);
            }
        }
    }
    return headToSeq;
}
