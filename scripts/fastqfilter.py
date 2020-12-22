#!/usr/bin/env python
"""
Author: Bruce Walker
"""


import sys
import argparse
import gzip
import bz2
import pysam
from Bio import SeqIO
import numpy as np

def openSeqFile(fileName):
    """
    Open a sequence file with SeqIO; can be fasta or fastq with optional gz or bz2 compression.
    Assumes fasta unless ".fastq" or ".fq" in the file name.
    :param fileName:
    :return: SeqIO.parse object
    """

    components = fileName.split('.')
    file = None
    if "gz" in components:
        file = gzip.open(fileName, 'rt')
    else:
        file = open(fileName, 'r')
    return SeqIO.parse(file, 'fastq')

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="output fastq file")
parser.add_argument("-b", "--bases", type=int, default=0, help="total bases desired")
parser.add_argument("-q", "--qual", type=float, default=0, help="minimum read quality")
parser.add_argument("-l", "--length", type=int, default=0, help="minimum read length")
parser.add_argument("-L", "--longlength", type=int, default=0, help="read length to include regardless")
parser.add_argument("-n", "--readnames", help="file containing read names to select")
parser.add_argument("input", help="input fastq file")
args = parser.parse_args()

total = 0
keepers = []
done = False
fastq = False

readnames = set()

if args.readnames:
    with open(args.readnames, 'r') as names:
        for name in names:
            readnames.add(name.strip())

for read in openSeqFile(args.input):
    if args.readnames and read.id not in readnames:
        continue

    length = len(read)
    if length < args.length:
        continue

    if done and length < args.longlength:
        continue

    if "phred_quality" in read.letter_annotations:
        fastq = True
        q = np.array(read.letter_annotations["phred_quality"])
        eq = -10.0 * np.log10(np.mean(np.power(10.0, -q/10.0)))
        if eq < args.qual:
            continue

    keepers.append(read)
    total += length

    done = args.bases and total >= args.bases
    
    if done and not args.longlength:
        break

if keepers:
    if args.output:
        output = open(args.output, 'w')
    else:
        output = sys.stdout

    SeqIO.write(keepers, output, "fastq" if fastq else "fasta")
