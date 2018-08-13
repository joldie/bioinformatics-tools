#!/usr/bin/python3

# number_fasta_sequences_v*.py
###################################
# Overview:
#   Simple script for adding numbers to sequence descriptions in FASTA files, in
#   order to make it easier to view in some analysis programs (e.g. UGENE)
#   Required input is any valid FASTA file with multiple sequences.
#
# Requirements:
#   biopython library: >> pip/pip3 install biopython
#
# License:
#   Public domain (Creative Commons Zero v1.0 Universal License)
#

# Imported libraries
######################
from Bio.Alphabet import IUPAC
from Bio import AlignIO
import os


# File data
######################
inputFileName = "file.fasta"
fileType = "fasta"
alphabetFormat = IUPAC.unambiguous_dna

inputFilePath = os.path.dirname(
    os.path.realpath(__file__)) + "/" + inputFileName
outputFilePath = inputFilePath + "_numbered." + fileType


######################

alignment = AlignIO.read(inputFilePath, fileType, alphabet=alphabetFormat)
numSequences = alignment.__len__()

for row in range(alignment.__len__()):
    alignment[row].description = "  __[" + \
        str(row+1) + "]__" + alignment[row].description

AlignIO.write(alignment, outputFilePath, fileType)
