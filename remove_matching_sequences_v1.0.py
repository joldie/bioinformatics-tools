#!/usr/bin/python3

# remove_matching_sequences_v*.py
###################################
# Overview:
#   Script analyses DNA multi-alignment sequence data and removes all identical
#   sequences found within.
#   User can adjust matching fraction (how "identical" the sequences must be).
#   Required input is any valid FASTA file with multiple sequences.
#   FASTA file outputted with "_cleaned" on end of filename.
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
from Bio.Align import MultipleSeqAlignment
import time
import os

# Measure execution time of script
start_time = time.time()


# Program settings
######################
matchingThreshold = 1  # 1 = identical (no mismatches)


# File data
######################
inputFileName = "file.fasta"
fileType = "fasta"
alphabetFormat = IUPAC.unambiguous_dna
gapCharacter = "-"

inputFilePath = os.path.dirname(
    os.path.realpath(__file__)) + "/" + inputFileName
outputFilePath = inputFilePath + "_cleaned." + fileType


######################

alignment = AlignIO.read(inputFilePath, fileType, alphabet=alphabetFormat)
sequenceLength = alignment.get_alignment_length()

if matchingThreshold < 0 or matchingThreshold > 1:
    print('Error: \'matchingThreshold\' must be between 0 and 1')
    exit()


def sequencesMatching(seq1, seq2):
    # Count number of columns where both sequences have valid data
    nonGapColumns = 0
    for column in range(sequenceLength):
        if seq1[column] != gapCharacter and seq2[column] != gapCharacter:
            nonGapColumns += 1

    # Determine maximum allowable mismatches (excluding gaps)
    maxAllowedMismatchedColumns = nonGapColumns - \
        int(nonGapColumns * matchingThreshold)

    # Count mismatched values between two sequences up until maximum reached
    mismatchedColumns = 0
    for column in range(sequenceLength):
        if seq1[column] != seq2[column]:
            mismatchedColumns += 1
        if mismatchedColumns > maxAllowedMismatchedColumns:
            return False
    return True


numberRemoved = 0
posSeq1 = 0
while posSeq1 < alignment.__len__():
    current_time = time.time() - start_time
    print("Start of sequence #" + str(posSeq1 + numberRemoved + 1) +
          " (current time = " + str(int(current_time)) + "s)")
    posSeq2 = posSeq1 + 1
    while posSeq2 < alignment.__len__():
        if sequencesMatching(alignment[posSeq1], alignment[posSeq2]):
            numberRemoved += 1
            print("Removing sequence #" + str(posSeq1 + numberRemoved +
                                              1) + " (" + alignment[posSeq2].name + ")")
            endSection = alignment[(posSeq2 + 1):]
            alignment = alignment[:(posSeq1 + 1)]
            alignment.extend(endSection)

        else:
            posSeq2 += 1
    posSeq1 += 1

print("Removed " + str(numberRemoved) +
      " sequences (" + str(alignment.__len__()) + " remaining)")

# Save "cleaned" file
AlignIO.write(alignment, outputFilePath, fileType)

# Print execution time of script
end_time = time.time()
total_time = end_time - start_time
print("Program execution time: " + str(int(total_time)) + "s")
