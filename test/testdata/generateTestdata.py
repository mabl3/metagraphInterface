import argparse
import json
import math
import os
import random
import re
import subprocess
import sys

parser = argparse.ArgumentParser(description = "Create testdata",
                                 formatter_class = argparse.RawTextHelpFormatter)
scriptArgs = parser.add_argument_group("Script Arguments")
scriptArgs.add_argument("--binsize",
                        dest = "binsize", 
                        metavar = "INT", 
                        type=int,
                        default=100,
                        help="binsize of metagraph")
scriptArgs.add_argument("--k",
                        dest = "k", 
                        metavar = "INT", 
                        type=int,
                        default=39,
                        help="k of metagraph")
scriptArgs.add_argument("--metagraph",
                        dest = "metagraph", 
                        metavar = "FILE", 
                        type=argparse.FileType("r"),
                        help="metagraph binary", 
                        required = True)
scriptArgs.add_argument("--nsequences",
                        dest = "nsequences", 
                        metavar = "INT", 
                        type=int,
                        default=10,
                        help="number of sequences to generate per species")
scriptArgs.add_argument("--nspecies",
                        dest = "nspecies", 
                        metavar = "INT", 
                        type=int,
                        default=10,
                        help="number of species to generate")
scriptArgs.add_argument("--sequence-lengths",
                        dest = "seqlen", 
                        metavar = "INT", 
                        type=int,
                        default=10000,
                        help="length of each sequence")
scriptArgs.add_argument("--similarity",
                        dest = "similarity", 
                        metavar = "FLOAT", 
                        type=float,
                        default=0.85,
                        help="sequence similarity between neighbouring species")
args = parser.parse_args()
metagraphBin = args.metagraph.name

# number of fasta files
numSpecies = args.nspecies
# number of sequences per fasta file
numSequencesPerSpecies = args.nsequences
# sequence similarity between sequence pairs
similarity = args.similarity
# lenths of the sequences
sequenceLength = args.seqlen
# k in graph
k = args.k
# bin size in graph
binsize = args.binsize

assert numSpecies >= 2
assert numSequencesPerSpecies >= 1
assert similarity > 0
assert similarity < 1
assert sequenceLength >= 1
assert k > 2
assert sequenceLength >= k
assert binsize >= 1


# find out which alphabet metagraph uses
alphabet = ""
caseSensitive = False
mg_build = os.path.dirname(metagraphBin)
mg_cc = os.path.join(mg_build, "compile_commands.json")
assert os.path.isfile(mg_cc), "Could not find metagraph compile_commands.json in " + mg_cc
with open(mg_cc, "r") as fh:
    metagraphCompileCommands = json.load(fh)
    for elem in metagraphCompileCommands:
        for key in elem:
            typestr = re.search("-D_[A-Z_]+_GRAPH", elem[key])
            if typestr:
                m = typestr.group(0)
                if m == "-D_DNA_GRAPH":
                    alphabet = "ACGT"
                elif m == "-D_DNA5_GRAPH":
                    alphabet = "ACGTN"
                elif m == "-D_DNA_CASE_SENSITIVE_GRAPH":
                    alphabet = "ACGTNacgt"
                    caseSensitive = True
                else:
                    sys.exit("Could not determine Metagraph alphabet")

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
        # if not DNA_CASE_SENSITIVE, omit softmask (as metagraph would convert lower case to upper case)
        insertBase = "N" if not caseSensitive and insertBase == newBase.lower() else insertBase

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
    sequences[0][sid] = baseSeq[0:-4] + "ACGT" # make sure seq ends in ACGT (not important here)
    for gid in range(1,numSpecies):
        sequences[gid][sid] = generateSimilarSequence(sequences[gid-1][sid], similarity)
        sequences[gid][sid] = sequences[gid][sid][0:-4] + "ACGT"

# make sure that at least one kmer occurs twice in the same bin
sequences[0][0] = sequences[0][0][0:k] + sequences[0][0][0:k] + sequences[0][0][(2*k):]
assert(len(sequences[0][0]) == sequenceLength)
# higher possibility that kmers are not present in each sequence
if k >= binsize and sequenceLength >= 3*binsize:
    sequences[1][0][2*binsize:3*binsize] = "N"*binsize

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
            # only add kmers that fit the alphabet
            if re.match("["+alphabet+"]{"+str(k)+"}", kmer):
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