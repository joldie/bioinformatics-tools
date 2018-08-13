#!/usr/bin/python3

# longest_identical_sequence_v*.py
###################################
# Overview:
#   Script analyses DNA multi sequence alignment data and searches for the
#   longest identical sequence within.
#   Required input is any valid FASTA file with multiple sequences.
#   CSV file outputted with detailed of longest identified sequences.
#
# Requirements:
#   biopython library: >> pip/pip3 install biopython
#
# License:
#   Public domain (Creative Commons Zero v1.0 Universal License)
#

# Measure execution time of script
import time
start_time = time.time()

# Input file details
directory = "/path/to/files"
filename = "file_name"

fasta_data = []
with open(directory + "/" + filename + '.fasta') as inputfile:
    for line in inputfile:
        fasta_data.append(line)

sequence_names = []
sequence_data = []
for row in fasta_data:
    if row[0] == ">":
        # Save sequence name, removing starting ">" and trailing newline "\n"
        sequence_names.append(row[1:len(row) - 1])
    else:
        # Save sequence data, removing trailing newline "\n"
        sequence_data.append(row[0:len(row) - 1])

# First line is reference sequence, to which all others will be compared
reference = sequence_data[0]

num_sequences = len(sequence_data)
num_genes = len(sequence_data[0])
min_compare_length = 2
result = []

for sequence in range(num_sequences):
    longest_subsequence = ""
    position_longest_subsequence = 0
    for test_length in range(min_compare_length, num_genes + 1):
        for position in range(num_genes - test_length + 1):
            compare1 = reference[position:position + test_length]
            compare2 = sequence_data[sequence][position:position + test_length]
            if "-" not in compare1 + compare2:
                if compare1 == compare2:
                    if test_length > len(longest_subsequence):
                        longest_subsequence = compare1
                        position_longest_subsequence = position
    if len(longest_subsequence) <= 10:
        snippet = longest_subsequence
    else:
        snippet = longest_subsequence[0:10] + "..."
    print(sequence_names[sequence] + ": Longest identical sequence = " +
          str(len(longest_subsequence)) + " bases")
    print("  [Position " + str(position_longest_subsequence) + "]: '" + snippet + "'")
    result.append([str(sequence), sequence_names[sequence], str(
        len(longest_subsequence)), str(position_longest_subsequence), snippet])

# Save results in CSV file
import csv
with open(directory + "/" + filename + '_identical_sequences.csv', 'w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=',')
    csv_writer.writerow(["Number", "Sequence name", "Length of longest identical sequence [bases]",
                         "Starting position of longest sequence", "Snippet of longest sequence"])
    for row in result:
        csv_writer.writerow(row)

# Print execution time of script
end_time = time.time()
total_time = end_time - start_time
print("Program execution time: " + str(int(total_time)) +
      "s (approx. " + str(round(total_time/num_sequences, 2)) + "s per sequence)")
