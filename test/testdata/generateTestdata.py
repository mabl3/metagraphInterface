import argparse
import json
import math
import random
import subprocess
import sys

parser = argparse.ArgumentParser(description = "Create testdata",
                                 formatter_class = argparse.RawTextHelpFormatter)
scriptArgs = parser.add_argument_group("Script Arguments")
scriptArgs.add_argument("--metagraph",
                        dest = "metagraph", 
                        metavar = "FILE", 
                        type=argparse.FileType("r"),
                        help="metagraph binary", 
                        required = True)
args = parser.parse_args()
metagraphBin = args.metagraph.name

# number of fasta files
numSpecies = 10
# number of sequences per fasta file
numSequencesPerSpecies = 10
# sequence similarity between sequence pairs
similarity = 0.85
# lenths of the sequences
sequenceLength = 10000
# k in graph
k = 39
# bin size in graph
binsize = 100

def generateSequence(seqLength):
    return("".join(random.choices(["A","C","G","T"], k=seqLength)))


def generateSimilarSequence(inputSeq, similarity):
    assert(similarity <= 1)
    assert(similarity > 0)

    seqList = []
    for c in inputSeq:
        p = [(1-similarity)/3]*4
        if c == "A":
            p[0] = similarity
        if c == "C":
            p[1] = similarity
        if c == "G":
            p[2] = similarity
        if c == "T":
            p[3] = similarity
        if c not in ["A","C","G","T"]:
            p = [0.25]*4

        newBase = random.choices(["A","C","G","T"], weights=p)[0]
        # from time to time, insert a softmask or unknown
        insertBase = random.choices([newBase, "N", newBase.lower()], weights=[0.99998, 0.00001, 0.00001])[0]
        seqList.append(insertBase)
    
    return("".join(seqList))

# generate sequences
sequences = dict()
for gid in range(numSpecies):
    sequences[gid] = dict()
    for sid in range(numSequencesPerSpecies):
        sequences[gid][sid] = ""

for sid in range(numSequencesPerSpecies):
    baseSeq = generateSequence(sequenceLength)
    sequences[0][sid] = baseSeq
    for gid in range(1,numSpecies):
        sequences[gid][sid] = generateSimilarSequence(sequences[gid-1][sid], similarity)

# make sure that at least one kmer occurs twice in the same bin
sequences[0][0] = sequences[0][0][0:k] + sequences[0][0][0:k] + sequences[0][0][(2*k):]
assert(len(sequences[0][0]) == sequenceLength)

# create data expected from graph
def genomeName(gid):
    return("genome"+str(gid)+".fa")

def sequenceName(sid):
    return("sequence"+str(sid))

kmers = dict()
for gid in sequences:
    for sid in sequences[gid]:
        seq = sequences[gid][sid]
        for i in range(0, len(seq)-k+1):
            kmer = seq[i:(i+k)]
            assert(len(kmer) == k)
            binIdx = math.floor(i/binsize) * binsize
            occ = tuple([genomeName(gid), sequenceName(sid), False, binIdx])
            if kmer not in kmers:
                kmers[kmer] = set() # set assumes that for kmers occuring >once in a bin, annotation is only reported once
                
            kmers[kmer].add(occ)

# need to convert sets and tuples to lists for json
for kmer in kmers:
    occset = kmers[kmer]
    kmers[kmer] = []
    for occ in occset:
        kmers[kmer].append(list(occ))

kmerNeighbours = dict()
for kmer in kmers:
    kmerNeighbours[kmer] = {"prev": [], "next": []}
    for b in ["A","C","G","T","N","a","c","g","t"]:
        prevKmer = b + kmer[0:(k-1)]
        nextKmer = kmer[1:k] + b
        assert(len(prevKmer) == k)
        assert(len(nextKmer) == k)
        if prevKmer in kmers:
            kmerNeighbours[kmer]["prev"].append(prevKmer)

        if nextKmer in kmers:
            kmerNeighbours[kmer]["next"].append(nextKmer)
    
    assert((len(kmerNeighbours[kmer]["prev"]) > 0) or (len(kmerNeighbours[kmer]["next"]) > 0))
    # also insert dummy nodes
    if len(kmerNeighbours[kmer]["prev"]) == 0:
        kmerNeighbours[kmer]["prev"] = ["$" + kmer[0:(k-1)]]
    
    if len(kmerNeighbours[kmer]["next"]) == 0:
        kmerNeighbours[kmer]["next"] = [kmer[1:k] + "$"]

# store data
with open("expectedKmers.json", "w") as fh:
    json.dump(kmers, fh)

with open("expectedKmerNeighbours.json", "w") as fh:
    json.dump(kmerNeighbours, fh)

# create the graph
for gid in sequences:
    with open(genomeName(gid), "w") as fh:
        for sid in sequences[gid]:
            fh.writelines(">"+sequenceName(sid)+"\n")
            fh.writelines(sequences[gid][sid]+"\n")
            fh.writelines("\n")
    
buildCommand = metagraphBin + " build -k " + str(k) + " --index-ranges 6 -o testdataGraph *.fa"
print("Running", buildCommand)
subprocess.run([buildCommand], shell = True, executable = '/bin/bash')

coordinateCommant = metagraphBin + " coordinate -i testdataGraph.dbg --coord-binsize " + str(binsize) + " *.fa"
print("Running", coordinateCommant)
subprocess.run([coordinateCommant], shell = True, executable = '/bin/bash')

sys.exit()


















































# generate C++ header with hardcoded data
content = []
content.append("#include <string>")
content.append("#include <unordered_map>")
content.append("")
content.append("#include <MetagraphInterface.h>")
content.append("")
content.append("")
content.append("")

### code chunks that initialize the NodeAnntoation vectors for each kmer
allKmerAnnotCode = []
for kmer in kmers:
    kmerAnnotCode = []
    kmerAnnotCode.append("    std::vector<MetagraphInterface::NodeAnnotation> annot"+kmer+"{")
    for annot in kmers[kmer]:
        kmerAnnotCode.append("      MetagraphInterface::NodeAnnotation{\""+annot[0]+"\", \""+annot[1]+"\", "+annot[2]+", "+str(annot[3])+"},")

    kmerAnnotCode[-1] = kmerAnnotCode[-1][0:(len(kmerAnnotCode[-1])-1)] # remove last komma
    kmerAnnotCode.append("    };")
    kmerAnnotCode.append("    expectedKmers.emplace(\""+kmer+"\", annot"+kmer+");")
    allKmerAnnotCode.append(kmerAnnotCode)

### create functions of 50 kmers each that fill an unordered_map of kmer to NodeAnnotation vector
for i in range(math.ceil(len(allKmerAnnotCode)/50)):
    content.append("void addExpectedKmers_"+str(i)+"(std::unordered_map<std::string, std::vector<MetagraphInterface::NodeAnnotation>> & expectedKmers) {")
    for c in range(i*50, min((i+1)*50, len(allKmerAnnotCode))):
        content.extend(allKmerAnnotCode[c])
        
    content.append("}")
    content.append("")

### 'parent' function to get an unordered_map of kmer to NodeAnnotation vector
content.append("auto expectedKmers() {")
content.append("    auto expectedKmers = std::unordered_map<std::string, std::vector<MetagraphInterface::NodeAnnotation>>{};")
for i in range(math.ceil(len(allKmerAnnotCode)/50)):
    content.append("    addExpectedKmers_"+str(i)+"(expectedKmers);")
content.append("}")
content.append("")
content.append("")
content.append("")

### code chunks that fill unordered_maps of kmers to previous or next kmer vectors
allKmerPredecessorCode = []
allKmerSuccessorCode = []
for kmer in kmerNeighbours:
    if len(kmerNeighbours[kmer]["prev"]) > 0:
        kmerPredecessorCode = ["    expectedPredecessors.emplace(\""+kmer+"\", std::vector<std::string>{"]
        for pred in kmerNeighbours[kmer]["prev"]:
            kmerPredecessorCode.append("      "+pred+",")
        kmerPredecessorCode[-1] = kmerPredecessorCode[-1][0:(len(kmerPredecessorCode[-1])-1)] # remove last komma
        kmerPredecessorCode.append("    };")
        allKmerPredecessorCode.append(kmerPredecessorCode)

    if len(kmerNeighbours[kmer]["next"]) > 0:
        kmerSuccessorCode = ["    expectedSuccessors.emplace(\""+kmer+"\", std::vector<std::string>{"]
        for succ in kmerNeighbours[kmer]["next"]:
            kmerSuccessorCode.append("      "+succ+",")
        kmerSuccessorCode[-1] = kmerSuccessorCode[-1][0:(len(kmerSuccessorCode[-1])-1)] # remove last komma
        kmerSuccessorCode.append("    };")
        allKmerSuccessorCode.append(kmerSuccessorCode)

### create functions of 50 kmers each that fill unordered_maps of kmer to kmer vectors
for i in range(math.ceil(len(allKmerPredecessorCode)/50)):
    content.append("void addPredecessorKmers_"+str(i)+"(std::unordered_map<std::string, std::vector<std::string>> & expectedPredecessors) {")
    for c in range(i*50, min((i+1)*50, len(allKmerPredecessorCode))):
        content.extend(allKmerPredecessorCode[c])
        
    content.append("}")
    content.append("")

content.append("auto expectedPredecessors() {")
content.append("    auto expectedPredecessors = std::unordered_map<std::string, std::vector<std::string>>{};")
for i in range(math.ceil(len(allKmerPredecessorCode)/50)):
    content.append("    addPredecessorKmers_"+str(i)+"(expectedPredecessors);")
content.append("}")
content.append("")
content.append("")
content.append("")

for i in range(math.ceil(len(allKmerSuccessorCode)/50)):
    content.append("void addSuccessorKmers_"+str(i)+"(std::unordered_map<std::string, std::vector<std::string>> & expectedSuccessors) {")
    for c in range(i*50, min((i+1)*50, len(allKmerSuccessorCode))):
        content.extend(allKmerSuccessorCode[c])
        
    content.append("}")
    content.append("")

content.append("auto expectedSuccessors() {")
content.append("    auto expectedSuccessors = std::unordered_map<std::string, std::vector<std::string>>{};")
for i in range(math.ceil(len(allKmerPredecessorCode)/50)):
    content.append("    addSuccessorKmers_"+str(i)+"(expectedSuccessors);")
content.append("}")
content.append("")

### write code to header file
for i in range(len(content)):
    content[i] = content[i] + "\n"
with open("loadExpectedTestdata.h", "w") as fh:
    fh.writelines(content)