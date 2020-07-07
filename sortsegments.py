#! /usr/bin/python3

"""
This script organizes and computes statistics about protein segments found
with findallsegments.py .

Author: Sébastien Légaré (ENS, INRIA, Paris, France)
Date: June 12, 2020
"""

# Minimum length of segments.
minlen = 6

# Output file.
outfile = open("sortedsegments.txt", "w")

# Input files.
segmentfiles = ["rawsegments/segments-r1ab_sars2.txt",
                "rawsegments/segments-spike_sars2.txt",
                #"rawsegments/segments-r1a_sars2.txt", 
                "rawsegments/segments-ns7a_sars2.txt",
                "rawsegments/segments-ap3a_sars2.txt",
                "rawsegments/segments-vme1_sars2.txt",
                "rawsegments/segments-ncap_sars2.txt",
                "rawsegments/segments-orf9b_sars2.txt",
                "rawsegments/segments-vemp_sars2.txt",
                "rawsegments/segments-ns6_sars2.txt",
                "rawsegments/segments-ns8_sars2.txt",
                "rawsegments/segments-ns7b_sars2.txt",
                "rawsegments/segments-y14_sars2.txt",
                "rawsegments/segments-a0a663dja2_sars2.txt"]

#  Viral protein descriptions.
viral_desc = {"P0DTD1": "R1AB_SARS2|YP_009724389.1 Replicase polyprotein 1ab",
              "P0DTC2": "SPIKE_SARS2|YP_009724390.1 Spike glycoprotein",
              "P0DTC1": "R1A_SARS2|YP_009725295.1 Replicase polyprotein 1a",
              "P0DTC7": "NS7A_SARS2|YP_009724395.1 Protein 7a",
              "P0DTC3": "AP3A_SARS2|YP_009724391.1 Protein 3a",
              "P0DTC5": "VME1_SARS2|YP_009724393.1 Membrane protein",
              "P0DTC9": "NCAP_SARS2|YP_009724397.2 Nucleoprotein",
              "P0DTD2": "ORF9B_SARS2|NONE Protein 9b",
              "P0DTC4": "VEMP_SARS2|YP_009724392.1 Envelope small membrane protein",
              "P0DTC6": "NS6_SARS2|YP_009724394.1 Non-structural protein 6",
              "P0DTC8": "NS8_SARS2|YP_009724396.1 Non-structural protein 8",
              "P0DTD8": "NS7B_SARS2|YP_009725318.1 Protein non-structural 7b",
              "P0DTD3": "Y14_SARS2|NONE Uncharacterized protein 14",
              "A0A663DJA2": "A0A663DJA2_SARS2|YP_009725255.1 ORF10 protein"}

# Read human protein descriptions.
fastafile = open("sequences/human_reviewed_canon.fasta", "r").readlines()
human_desc = {}
for line in fastafile:
    if ">sp" in line:
        pipe1 = line.index("|")
        pipe2 = line.rfind("|")
        org = line.rfind("OS=")
        uniprotac = line[pipe1+1:pipe2]
        desc = line[pipe2+1:org-1]
        human_desc[uniprotac] = desc

# Gather sequence segments.
segments = {}
for segmentfile in segmentfiles:
    segcontent = open(segmentfile, "r").readlines()
    for line in segcontent:
        tokens = line.split()
        sequence = tokens[0]
        if len(sequence) >= minlen:
            pipe1 = tokens[1].index("|")+1
            pipe2 = tokens[1][pipe1:].index("|")
            id1 = tokens[1][pipe1:pipe1+pipe2]
            start1 = int(tokens[2][1:-1])+1
            end1 = int(tokens[3][:-1])+1
            pipe3 = tokens[4].index("|")+1
            pipe4 = tokens[4][pipe3:].index("|")
            id2 = tokens[4][pipe3:pipe3+pipe4]
            start2 = int(tokens[5][1:-1])+1
            end2 = int(tokens[6][:-1])+1
            range1 = "{}-{}".format(start1, end1)
            range2 = "{}-{}".format(start2, end2)
            if sequence not in segments.keys():
                segments[sequence] = {"viral": {id1: [range1]},
                                      "human": {id2: [range2]}}
            else:
                if id1 not in segments[sequence]["viral"].keys():
                    segments[sequence]["viral"][id1] = [range1]
                elif range1 not in segments[sequence]["viral"][id1]:
                    segments[sequence]["viral"][id1].append(range1)
                if id2 not in segments[sequence]["human"].keys():
                    segments[sequence]["human"][id2] = [range2]
                elif range2 not in segments[sequence]["human"][id2]:
                    segments[sequence]["human"][id2].append(range2)

# Count number of different regions on proteins.
for sequence in segments.keys():
    nvir = 0
    for id1 in segments[sequence]["viral"].keys():
        nvir += len(segments[sequence]["viral"][id1])
    segments[sequence]["viral"]["n"] = nvir
    nhum = 0
    for id2 in segments[sequence]["human"].keys():
        nhum += len(segments[sequence]["human"][id2])
    segments[sequence]["human"]["n"] = nhum
    segments[sequence]["n"] = nvir + nhum

# Transform segment dict into a list to allow sorting.
seglist = []
for sequence in segments.keys():
    new_seg = {"seq": sequence,
               "viral": segments[sequence]["viral"],
               "human": segments[sequence]["human"],
               "n": segments[sequence]["n"],
               "l": len(sequence)}
    seglist.append(new_seg)

# Sort by segment length and number of times found in proteins.
sorted_segments1 = sorted(seglist, key=lambda x: x["n"], reverse=True)
sorted_segments = sorted(sorted_segments1, key=lambda x: x["l"], reverse=True)

# Write resulting segments.
for seg in sorted_segments:
    outfile.write("Sequence {},  ".format(seg["seq"]))
    outfile.write("Len: {}   Num: {} ".format(seg["l"], seg["n"]))
    outfile.write("({} viral, ".format(seg["viral"]["n"]))
    outfile.write("{} human)\n".format(seg["human"]["n"]))
    outfile.write("Viral: (")
    outfile.write("{} segments among ".format(seg["viral"]["n"]))
    outfile.write("{} proteins)\n".format(len(seg["viral"].keys())-1))
    for ptn in seg["viral"].keys():
        if ptn != "n":
            outfile.write("{} ".format(ptn))
            outfile.write("{}  ".format(viral_desc[ptn]))
            outfile.write("{}\n".format(" ".join(seg["viral"][ptn])))
    outfile.write("Human: (")
    outfile.write("{} segments among ".format(seg["human"]["n"]))
    outfile.write("{} proteins)\n".format(len(seg["human"].keys())-1))
    for ptn in seg["human"].keys():
        if ptn != "n":
            outfile.write("{} ".format(ptn))
            outfile.write("{}  ".format(human_desc[ptn]))
            outfile.write("{}\n".format(" ".join(seg["human"][ptn])))
    outfile.write("----------\n")

