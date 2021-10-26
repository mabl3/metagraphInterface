"""
    Create a metagraph file, an annotation file and a meta file.

    The graph and annotation files can be loaded in C++ with `MetagraphInterface`,
    each k-mer is annotated with its sequence headers and the exact positions.

    The meta file stores additional information that is not (completely) encoded
    in the graph. It contains a mapping of the filenames to the sequence headers 
    that are present in that respective file, a mapping of each sequence header 
    to the respective filename, and a mapping of sequence headers to sequence
    lengths. The mappings can be used in seedFinding to quickly prepare a run. The
    file is in JSON format.

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
scriptArgs.add_argument("--only-meta",
                        dest = "onlyMeta",
                        action = 'store_true',
                        help = "Do not create graph, only meta file")
scriptArgs.add_argument("--output",
                        dest = "basename", 
                        metavar = "Metagraph Basefilename", 
                        type = str,
                        help = "Basename to use for metagraph files", 
                        required = True)
scriptArgs.add_argument("--tmpdir",
                        dest = "tmpdir", 
                        metavar = "Path", 
                        type = str,
                        default = "/tmp/",
                        help = "Temporary directory for intermediate files")
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

if not args.onlyMeta:
    # create temporary working directory
    tmpdir = os.path.join(args.tmpdir, "createGraph_" + secrets.token_hex(7))
    tmpAnnoDir = os.path.join(tmpdir, "annotations")
    assert not os.path.exists(tmpdir), tmpdir + " already exists"
    os.mkdir(tmpdir)
    os.mkdir(tmpAnnoDir)
    assert os.path.exists(tmpdir), tmpdir + " could not be created"
    assert os.path.exists(tmpAnnoDir), tmpAnnoDir + " could not be created"

    # create graph
    tmpfilebase = os.path.join(tmpdir, outfilebase)
    buildCommand = metagraphBin + " build -k " + str(args.k) + " --index-ranges 6 -p " + str(multiprocessing.cpu_count()) + " -o " + tmpfilebase + " " + (" ".join(fastas))
    print("[INFO] >>> Running", buildCommand)
    subprocess.run([buildCommand], shell = True, executable = '/bin/bash')

# create annotations
seenSeqs = set() # check if all sequence headers are unique, warn if not
seqToFile = {}
fileToSeq = {}
seqToLen = {}
i = 0
for fasta in fastas:
    _, fastaname = os.path.split(fasta)
    genome, _ = os.path.splitext(fastaname)
    if genome not in fileToSeq:
        fileToSeq[genome] = []

    for seq in SeqIO.parse(fasta, "fasta"):
        if seq.id in seenSeqs:
            print("[WARNING] >>> Sequence ID '"+seq.id+"' is ambiguous!")

        else:
            if not args.onlyMeta:
                tmpfile = os.path.join(tmpAnnoDir, str(i)+".fa")
                with open(tmpfile, 'wt') as fh:
                    SeqIO.write(seq, fh, 'fasta')

                annoCommand = metagraphBin + " annotate --coordinates --anno-header -i " + tmpfilebase + ".dbg -o " + os.path.join(tmpAnnoDir, outfilebase + "_"+str(i)) + " " + tmpfile
                print("[INFO] >>> Running", annoCommand)
                subprocess.run([annoCommand], shell = True, executable = '/bin/bash')

            seenSeqs.add(seq.id)
            fileToSeq[genome].append(seq.id)
            seqToFile[seq.id] = genome
            seqToLen[seq.id] = len(seq)
            i += 1


if not args.onlyMeta:
    transformCommand = metagraphBin + " transform_anno --anno-type column_coord -p " + str(multiprocessing.cpu_count()) + " -o " + tmpfilebase + " " + os.path.join(tmpAnnoDir, "*.column.annodbg")
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

with open(args.basename+".meta.json", 'wt') as fh:
    json.dump({'fileToSeq': fileToSeq,
               'seqToFile': seqToFile,
               'seqToLen': seqToLen}, fh)

print("[INFO] >>> Done.")