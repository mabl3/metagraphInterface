"""
    Create a metagraph file, an annotation file, a sequence header mapping and position correction.

    The graph and annotation files can be loaded in C++ with `MetagraphInterface`,
    each k-mer is annotated with its sequence headers and the exact positions.

    The mappings and position correction are stored in two JSON files, 
    one maps the filename to the sequence headers that are present in that file, 
    the other maps each sequence header to the respective filename and the correction
    maps each sequence header to a position correction value, needed for MetagraphInterface. 
    The mappings can be used in seedFinding to prepare a valid IdentifierMapping.

    If sequence headers are ambiguous, the mapping will be incomplete and sequences
    will not be represented correctly!
"""

import argparse
from Bio import SeqIO
import json
import multiprocessing
import os
import secrets
import shutil
import subprocess

parser = argparse.ArgumentParser(description = "Create Metagraph",
                                 formatter_class = argparse.RawTextHelpFormatter)
scriptArgs = parser.add_argument_group("Script Arguments")
scriptArgs.add_argument("--input",
                        dest = "input", 
                        metavar = "LIST OF FILES", 
                        type = argparse.FileType("r"), 
                        nargs = "+",
                        help = "Filenames of the input fasta files",
                        required = True)
scriptArgs.add_argument("--k",
                        dest = "k", 
                        metavar = "INT", 
                        type = int,
                        help = "k of metagraph",
                        required = True)
scriptArgs.add_argument("--metagraph",
                        dest = "metagraph", 
                        metavar = "FILE", 
                        type = argparse.FileType("r"),
                        help = "metagraph binary", 
                        required = True)
scriptArgs.add_argument("--output",
                        dest = "basename", 
                        metavar = "Metagraph Basefilename", 
                        type = str,
                        help = "Basename to use for metagraph files", 
                        required = True)
args = parser.parse_args()
metagraphBin = args.metagraph.name
fastas = [f.name for f in args.input]

print("[INFO] >>> Preparing Metagraph Build")

# check output dir
outpath, outfilebase = os.path.split(args.basename)
if outpath == '':
    outpath = '.'

outpath = os.path.abspath(outpath)
assert os.path.isdir(outpath), outpath + " is not a directory" 

# create temporary working directory
tmpdir = "/tmp/createGraph_" + secrets.token_hex(7)
assert not os.path.exists(tmpdir), tmpdir + " already exists"
os.mkdir(tmpdir)

# create graph
seenSeqs = set() # check if all sequence headers are unique, warn if not
seqToFile = {}
fileToSeq = {}
seqToCorrection = {}
posCorrection = 0
for fasta in fastas:
    _, fastaname = os.path.split(fasta)
    genome, _ = os.path.splitext(fastaname)
    if genome not in fileToSeq:
        fileToSeq[genome] = []

    for seq in SeqIO.parse(fasta, "fasta"):
        if seq.id in seenSeqs:
            print("[WARNING] >>> Sequence ID '"+seq.id+"' is ambiguous!")

        else:
            seenSeqs.add(seq.id)
            fileToSeq[genome].append(seq.id)
            seqToFile[seq.id] = genome
            seqToCorrection[seq.id] = posCorrection
            posCorrection += (len(seq) - int(args.k) + 1)



tmpfilebase = os.path.join(tmpdir, outfilebase)
buildCommand = metagraphBin + " build -k " + str(args.k) + " --index-ranges 6 -p " + str(multiprocessing.cpu_count()) + " -o " + tmpfilebase + " " + (" ".join(fastas))
print("[INFO] >>> Running", buildCommand)
subprocess.run([buildCommand], shell = True, executable = '/bin/bash')

annoCommand = metagraphBin + " annotate --coordinates --anno-header --separately -p " + str(multiprocessing.cpu_count()) + " -i " + tmpfilebase + ".dbg -o " + tmpfilebase + "Annotations " + (" ".join(fastas))
print("[INFO] >>> Running", annoCommand)
subprocess.run([annoCommand], shell = True, executable = '/bin/bash')

transformCommand = metagraphBin + " transform_anno --anno-type column_coord -p " + str(multiprocessing.cpu_count()) + " -o " + tmpfilebase + " " + tmpfilebase + "Annotations/*.column.annodbg"
print("[INFO] >>> Running", transformCommand)
subprocess.run([transformCommand], shell = True, executable = '/bin/bash')

assert os.path.isfile(tmpfilebase+".dbg"), tmpfilebase+".dbg does not exist"
assert os.path.isfile(tmpfilebase+".column_coord.annodbg"), tmpfilebase+".column_coord.annodbg does not exist"

# move graph to destination
print("[INFO] >>> Finishing Build")
shutil.move(tmpfilebase+".dbg", outpath)
shutil.move(tmpfilebase+".column_coord.annodbg", outpath)

assert os.path.isfile(args.basename+".dbg"), args.basename+".dbg does not exist"
assert os.path.isfile(args.basename+".column_coord.annodbg"), args.basename+".column_coord.annodbg does not exist"

shutil.rmtree(tmpdir, onerror=lambda:print("Failed to remove", tmpdir))

with open(args.basename+".seqToFile.json", 'wt') as fh:
    json.dump(seqToFile, fh)

with open(args.basename+".fileToSeq.json", 'wt') as fh:
    json.dump(fileToSeq, fh)

with open(args.basename+".positionCorrection.txt", 'wt') as fh:
    fh.writelines([seq+"\n"+str(seqToCorrection[seq])+"\n\n" for seq in seqToCorrection])

print("[INFO] >>> Done.")