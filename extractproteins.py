#! /usr/bin/python3

"""
This script organizes and computes statistics about protein segments found
with findallsegments.py .

Author: Sébastien Légaré (ENS, INRIA, Paris, France)
Date: June 12, 2020
"""

# Minimum length of segments.
minlen = 6

# List of proteins to include in extraction.
# If selectionfile in None, keep all proteins.
selectionfile = None
#selectionfile = "proteinlists/uniprotac-endothelial.txt"


# Output file.
outfile = open("proteinmatches.txt", "w")
#outfile = open("proteinmatches-endothelial.txt", "w")


# Input file.
segmentfile = open("sortedsegments.txt", "r").readlines()

# Read selection of proteins.
ptnsel = []
if selectionfile != None:
    selectioncontent = open(selectionfile, "r").readlines()
    for line in selectioncontent:
        if len(line) > 3:
            ptnsel.append(line[:-1])

# Gather proteins.
humanproteins = {}
virallist = []
viral_desc = {}
human_desc = {}
gather = None
for i in range(len(segmentfile)):
    line = segmentfile[i]
    tokens = line.split()
    if "Sequence" in line:
        seq = tokens[1][:-1]
    if gather == "viral" and "Human:" not in line:
        first_space = line.index(" ")
        double_space = line.index("  ")
        viral_ptn = line[:first_space]
        ranges = line[double_space+2:-1]
        new_viral = {"ptn": viral_ptn, "ranges": ranges}
        virallist.append(new_viral)
        desc = line[first_space+1:double_space]
        if viral_ptn not in viral_desc.keys():
            viral_desc[viral_ptn] = desc
    if gather == "human" and "----------" not in line:
        first_space = line.index(" ")
        double_space = line.index("  ")
        human_ptn = line[:first_space]
        ranges = line[double_space+2:-1]
        new_seg = {"seq": seq, "ranges": ranges, "viral": virallist}
        if len(ptnsel) == 0 or (len(ptnsel) > 0 and human_ptn in ptnsel):
            if len(seq) >= minlen:
                if human_ptn not in humanproteins.keys():
                    humanproteins[human_ptn] = [new_seg]
                else:
                    humanproteins[human_ptn].append(new_seg)
        desc = line[first_space+1:double_space]
        if human_ptn not in human_desc.keys():
            human_desc[human_ptn] = desc
    if "Viral:" in line:
        gather = "viral"
    if "Human:" in line:
        gather = "human"
    if "-----" in line:
        gather = None
        virallist = []

# Order segments within human proteins.
proteins_with_sorted_segs = {}
for ptn in humanproteins.keys():
    seglist = humanproteins[ptn]
    for seg in seglist:
        range_tokens = seg["ranges"].split()
        first_range = range_tokens[0]
        dash = first_range.index("-")
        first_residue = int(first_range[:dash])
        seg["first"] = first_residue
    sorted_segments = sorted(seglist, key=lambda x: x["first"])
    proteins_with_sorted_segs[ptn] = sorted_segments

# Transform protein dict into a list to allow sorting.
proteinlist = []
for ptn in proteins_with_sorted_segs.keys():
    nranges = 0
    seq_lens = []
    for seg in proteins_with_sorted_segs[ptn]:
        range_tokens = seg["ranges"].split()
        nranges += len(range_tokens)
        seq_lens.append(len(seg["seq"]))
    new_ptn = {"ptn": ptn,
               "segments": proteins_with_sorted_segs[ptn],
               "nranges": nranges,
               "longest": max(seq_lens)}
    proteinlist.append(new_ptn)

# Sort proteins by highest number of segments.
sortedproteins1 = sorted(proteinlist, key=lambda x: x["nranges"], reverse=True)
sortedproteins = sorted(sortedproteins1, key=lambda x: x["longest"], reverse=True)

# Write resulting segments.
for ptn in sortedproteins:
    outfile.write(">{} {} ".format(ptn["ptn"], human_desc[ptn["ptn"]]))
    outfile.write("(N segments = {})\n".format(ptn["nranges"]))
    for seg in ptn["segments"]:
        outfile.write("{} {}   Viral matches: ".format(seg["seq"],
                                                       seg["ranges"]))
        for i in range(len(seg["viral"])):
            viral = seg["viral"][i]
            outfile.write("{} {} ".format(viral["ptn"],
                                         viral_desc[viral["ptn"]]))
            outfile.write("{}".format(viral["ranges"]))
            if i < len(seg["viral"])-1:
                outfile.write(" / ")
        outfile.write("\n")
    outfile.write("----------\n")

#for p in sortedproteins:
#    print(p)

