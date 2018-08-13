#!/usr/bin/python3

# de-number_fasta_sequences_v*.py
###################################
# Overview:
#   Simple script for removing numbers to sequence descriptions in FASTA files.
#   Effectively undoes number_fasta_sequences_v*.py script.
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
outputFilePath = inputFilePath + "_de-numbered." + fileType


######################

alignment = AlignIO.read(inputFilePath, fileType, alphabet=alphabetFormat)
numSequences = alignment.__len__()

for row in range(alignment.__len__()):
    stop = alignment[row].description.index("]__") + 3
    alignment[row].description = alignment[row].description[stop:]

AlignIO.write(alignment, outputFilePath, fileType)
