#! /usr/bin/python3

"""
This script finds common segments between two sequences, no matter the order in
which they are found in the complete sequence.

Author: Sébastien Légaré (ENS, INRIA, Paris, France)
Date: June 12, 2020
"""

import time

## Input files.
## (Each run in separate execution because some take several hours)
fastafile = "sequences/r1ab_sars2.fasta"       #  1
#fastafile = "sequences/spike_sars2.fasta"      #  2
#fastafile = "sequences/r1a_sars2.fasta"        #  3
#fastafile = "sequences/ns7a_sars2.fasta"       #  4
#fastafile = "sequences/ap3a_sars2.fasta"       #  5
#fastafile = "sequences/vme1_sars2.fasta"       #  6
#fastafile = "sequences/ncap_sars2.fasta"       #  7
#fastafile = "sequences/orf9b_sars2.fasta"      #  8
#fastafile = "sequences/vemp_sars2.fasta"       #  9
#fastafile = "sequences/ns6_sars2.fasta"        # 10
#fastafile = "sequences/ns8_sars2.fasta"        # 11
#fastafile = "sequences/ns7b_sars2.fasta"       # 12
#fastafile = "sequences/y14_sars2.fasta"        # 13
#fastafile = "sequences/a0a663dja2_sars2.fasta" # 14

fastafile2 = "sequences/human_reviewed_canon.fasta"

# Minimum length of segments.
minlen = 5


# Function definitions.
def read_fasta_file(fasta):
    """ Get the content of a fasta file. """

    ptn_list = []
    fasta_content = open(fasta, "r")
    new_ptn = None
    for line in fasta_content:
        if ">sp" in line or ">tr" in line:
            if new_ptn != None:
                new_ptn["seq"] = sequence
                ptn_list.append(new_ptn)
            tokens = line.split()
            new_ptn = {"id": tokens[0] }
            sequence = ""
        else:
            sequence += line[:-1]
    new_ptn["seq"] = sequence
    ptn_list.append(new_ptn)

    return ptn_list


def findsegments(id1, seq1, id2, seq2, minlen):
    """
    Find common segments between two sequences, no matter the
    order in which they are found in the complete sequence.
    """

    segments = ""

    # Initialize list of corresponding residues.
    correspondances = []
    for res in seq1:
        correspondances.append([])
    
    # Main loop.
    for i in range(len(seq1)-minlen):
        seg1 = seq1[i:i+minlen]
        for j in range(len(seq2)-minlen):
            if j not in correspondances[i]:
                seg2 = seq2[j:j+minlen]
                if seg1 == seg2:
                    # Look if the segment is longer than minlen.
                    segments_equal = True
                    prev1 = seg1
                    prev2 = seg2
                    extend = 1
                    while segments_equal == True:
                        i_end = i+minlen+extend
                        j_end = j+minlen+extend
                        ext1 = seq1[i:i_end]
                        ext2 = seq2[j:j_end]
                        if i_end > len(seq1) or j_end > len(seq2):
                            seqend = True
                        else:
                            seqend = False
                        if ext1 != ext2 or seqend == True:
                            segments_equal = False
                            segments += "{}   ".format(prev1)
                            segments += "{} [{}, {}] ".format(id1, i, i_end-2)
                            segments += "  "
                            segments += "{} [{}, {}] ".format(id2, j, j_end-2)
                            segments += "\n"
                            # Add residues to correspondance list.
                            for k in range(minlen+extend-1):
                                l = i+k
                                m = j+k
                                correspondances[l].append(m)
                        prev1 = ext1
                        prev2 = ext2
                        extend += 1

    return segments


# Run code.

# Open output file.
slash = fastafile.index("/")
dot = fastafile.index(".")
viral_seq = fastafile[slash+1:dot]
outfile = "rawsegments/segments-{}.txt".format(viral_seq)
outputfile = open(outfile, "w")

# Read fasta files.
ptn_list1 = read_fasta_file(fastafile)
ptn_list2 = read_fasta_file(fastafile2)

time_start = time.perf_counter()

# Find segments.
counter = 1
for ptn1 in ptn_list1:
    print(ptn1["id"])
    for ptn2 in ptn_list2:
        segments_str = findsegments(ptn1["id"], ptn1["seq"],
                                    ptn2["id"], ptn2["seq"], minlen)
        outputfile.write(segments_str)
        counter += 1
        if counter%100 == 0:
            print("{} / {}".format(counter, len(ptn_list2)))

time_stop = time.perf_counter()
time_diff = time_stop - time_start
print("\nCalculation time : {:.2f}s\n".format(time_diff))

