#!/usr/bin/python3

# find_conserved_blocks_v*.py
###########################
# Overview:
#   Script analyses DNA multi sequence alignment data and searches for areas
#   (blocks) which are conserved (i.e. identical or very similar). Required
#   input is any valid FASTA file with multiple sequences.
#   User can adjust matching fraction (1 = identical match) and length of
#   matching region. Computational progress shown in console and text file
#   generated at end showing results.
#
# Requirements:
#   biopython library: >> pip/pip3 install biopython
#
# License:
#   Public domain (Creative Commons Zero v1.0 Universal License)
#

# Imported libraries
###########################
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import time
import os

# Measure execution time of script
start_time = time.time()


# Program settings
###########################
matchingFraction = 0.85  # 1 = identical (no mismatches)
minLength = 200          # minimum number of base pairs matching in sequence
minSequences = 7        # minimum number of sequences in a block with matching sequence


# Input file data
###########################
inputFileName = "file.fasta"
inputFileType = "fasta"
alphabetFormat = IUPAC.unambiguous_dna
gapCharacter = "-"

inputFilePath = os.path.dirname(
    os.path.realpath(__file__)) + "/" + inputFileName
outputFilePath = inputFilePath + "_(" + str(matchingFraction) + "-" + \
    str(minLength) + "-" + str(minSequences) + ").txt"


###########################

alignment = AlignIO.read(inputFilePath, inputFileType, alphabet=alphabetFormat)
sequenceLength = alignment.get_alignment_length()
numSequences = alignment.__len__()

# Error checking for user input
if matchingFraction < 0 or matchingFraction > 1:
    print('Error: \'matchingFraction\' must be between 0 and 1')
    exit()

if minLength >= sequenceLength or minLength <= 0:
    print('Error: \'minLength\' must be a positive number less than the sequence length in input file')
    exit()

if minSequences >= numSequences or minSequences <= 0:
    print('Error: \'minSequences\' must be a positive number less than the number of sequences in input file')
    exit()


def isBlockConserved(blockSequence):
    maxAllowedMismatchedColumns = blockSequence.get_alignment_length(
    ) - int(blockSequence.get_alignment_length() * matchingFraction)
    mismatchedColumns = 0
    for column in range(blockSequence.get_alignment_length()):
        # First, check if entire first column is filled with gap characters
        # If so, return False
        if column == 0:
            gapCharacterCount = 0
            for row in range(blockSequence.__len__()):
                if blockSequence[row][column] == gapCharacter:
                    gapCharacterCount += 1
            if gapCharacterCount == blockSequence.__len__():
                return False
        # Find the first non-gap character in column
        firstNonGapCharacter = ""
        for row in range(blockSequence.__len__()):
            if blockSequence[row][column] != gapCharacter and firstNonGapCharacter == "":
                firstNonGapCharacter = blockSequence[row][column]
                break
        # Next, count number of columns with non-identical values or gap characters
        # If count increases above threshold, return False
        for row in range(blockSequence.__len__()):
            if blockSequence[row][column] == gapCharacter or blockSequence[row][column] != firstNonGapCharacter:
                mismatchedColumns += 1
                if mismatchedColumns > maxAllowedMismatchedColumns:
                    return False
                break
    # If loop executes completely to end, return True
    return True


def isBlockAlreadyCounted(startRow, startColumn, endRow, endColumn):
    for block in conservedBlocksInfo:
        if startRow >= block[0] and startColumn >= block[1] and endRow <= block[2] and endColumn <= block[3]:
            return True
    return False


def removeTrailingBlankColumns(blockSequence):
    removedColumns = 0
    blankColumn = True
    while blankColumn:
        gapCharacterCount = 0
        for row in range(blockSequence.__len__()):
            if blockSequence[row][-1] == gapCharacter:
                gapCharacterCount += 1
        if gapCharacterCount == blockSequence.__len__():
            # Remove last column, as it contains only gap characters
            blockSequence = blockSequence[:, :-1]
            removedColumns += 1
        else:
            blankColumn = False
    return blockSequence, removedColumns


def getNumberBlankColumns(blockSequence):
    numberColumns = 0
    for column in range(blockSequence.get_alignment_length()):
        gapCharacterCount = 0
        for row in range(blockSequence.__len__()):
            if blockSequence[row][column] == gapCharacter:
                gapCharacterCount += 1
        if gapCharacterCount == blockSequence.__len__():
            numberColumns += 1
    return numberColumns


def getMismatchString(blockSequence):
    mismatchString = ""
    for column in range(blockSequence.get_alignment_length()):
        mismatch = False
        firstNonGapCharacter = ""
        for row in range(blockSequence.__len__()):
            if blockSequence[row][column] != gapCharacter and firstNonGapCharacter == "":
                firstNonGapCharacter = blockSequence[row][column]
                break
        for row in range(blockSequence.__len__()):
            if blockSequence[row][column] != firstNonGapCharacter:
                mismatchString += "X"
                mismatch = True
                break
        if not mismatch:
            mismatchString += "-"
    return mismatchString


def getNumberPositionString(length, start, multiple):
    returnString = ""
    i = start
    while i < length + start:
        if i % multiple == 0:
            returnString += str(i)
            # Move counter forward by number of extra characters used to print string
            i += (len(str(i)) - 1)
        else:
            returnString += " "
        i += 1
    return returnString


print("Input file = " + inputFilePath)
print("\'matchingFraction\' = " + str(matchingFraction))
print("\'minLength\' = " + str(minLength))
print("\'minSequences\' = " + str(minSequences))
print("")

conservedBlocks = []
conservedBlocksInfo = []

for startSequence in range(numSequences - minSequences + 1):
    current_time = time.time() - start_time
    print("Start of sequence #" + str(startSequence + 1) +
          " (current time = " + str(int(current_time)) + "s)")
    for index in range(sequenceLength - minLength + 1):
        length = minLength
        sequences = minSequences
        # Only test block if it (or larger version of it) not already added
        if not isBlockAlreadyCounted(startSequence, index, startSequence + sequences, index + length):
            largestConservedBlock = MultipleSeqAlignment([])  # Empty object

            consensusRight = True
            consensusDown = True

            while consensusRight:
                block = alignment[startSequence:(
                    startSequence + sequences), index:(index + length)]
                if (index + length) > sequenceLength:
                    # Reached end of sequence
                    consensusRight = False
                elif not isBlockConserved(block):
                    consensusRight = False
                    if length == minLength:
                        consensusDown = False  # Don't bother searching down, as start block is not conserved
                    else:
                        length -= 1
                else:
                    length += 1
                    largestConservedBlock = block

            while consensusDown:
                block = alignment[startSequence:(
                    startSequence + sequences), index:(index + length)]
                if (startSequence + sequences) > numSequences:
                    # Reached end of sequence list
                    consensusDown = False
                elif not isBlockConserved(block):
                    consensusDown = False
                    sequences -= 1
                else:
                    sequences += 1
                    largestConservedBlock = block

            if index != sequenceLength and length > minLength:
                # Before saving data, remove any columns of gap characters at end of block
                largestConservedBlock, removedColumns = removeTrailingBlankColumns(
                    largestConservedBlock)
                # Save block data in list
                conservedBlocks.append(largestConservedBlock)
                conservedBlocksInfo.append(
                    [startSequence, index, startSequence + sequences, index + length])
                # Print information about saved data
                print("Block found: " + str(conservedBlocks[-1].__len__()) + " sequences, " + str(
                    conservedBlocks[-1].get_alignment_length()) + " long (from position " + str(index + 1) + " to " + str(index + length + 1) + ")")


# Print output to text file
textFile = open(outputFilePath, "w")
textFile.write("Python script = " + os.path.basename(__file__) + "\n")
textFile.write("Input file = " + inputFilePath + "\n")
textFile.write("\'matchingFraction\' = " + str(matchingFraction) + "\n")
textFile.write("\'minLength\' = " + str(minLength) + "\n")
textFile.write("\'minSequences\' = " + str(minSequences) + "\n\n")

i = 0
for block in conservedBlocks:
    textFile.write("Number blank columns = " +
                   str(getNumberBlankColumns(block)) + "\n")
    textFile.write("Length = " + str(block.get_alignment_length() - getNumberBlankColumns(block)) + " / Length with gaps = " + str(block.get_alignment_length()
                                                                                                                                   ) + " (position " + str(conservedBlocksInfo[i][1] + 1) + " to " + str(conservedBlocksInfo[i][1] + block.get_alignment_length() + 1) + ")\n")
    textFile.write("Number of sequences = " + str(block.__len__()) + "\n")
    for row in range(block.__len__()):
        textFile.write("[" + str(conservedBlocksInfo[i][0] +
                                 row + 1) + "] " + block[row].name + "\n")
    textFile.write(getNumberPositionString(
        block.get_alignment_length(), 1, 50) + "\n")
    textFile.write(getMismatchString(block) + "\n")
    for row in range(block.__len__()):
        textFile.write(str(block[row, :].seq) + "\n")
    textFile.write("\n")
    i += 1

textFile.close()


# Print execution time of script
end_time = time.time()
total_time = end_time - start_time
print("Program execution time: " + str(int(total_time)) + "s")
