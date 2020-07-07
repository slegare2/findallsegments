# findallsegments

This code finds identical continuous segments between two proteins. It is here used to compare SARS-CoV-2 and human proteins.

## Usage

The main script is `findallsegments.py`. It reads the viral and human sequences from directory `sequences` and writes
the identical segments found for each viral protein in directory `rawsegments`.

Function `findsegments` in `findallsegments.py` is where the systematic window sliding search is done.

The two remaining scripts, `sortsegments.py`  and `extractproteins.py`, are mostly about organizing the results for better readability.

Script `sortsegments.py` sorts the segments obtained by length and occurrence and writes the result to file `sortedsegments.txt`.

The script `extractproteins.py` rewrites the results to file `proteinmatches.txt` in a protein-centric fashion rather that segment-centric.
This last script accepts an optional list of proteins. If a list is provided, proteins that are not in the list are excluded from
resulting matches.

## Details on the protein sequences used

Human sequences file `sequences/human_reviewed_canon.fasta` was downloaded from UniProt on June 3, 2020
using the following query: `organism:"Homo sapiens (Human) [9606]" AND reviewed:yes`.

SARS-Cov-2 viral sequences files `sequences/r1ab_sars2.fasta`, `sequences/spike_sars2.fasta`, etc,
were downloaded from UniProt on June 3, 2020 using the following query: `organism:sars-cov-2 AND reviewed:yes`.

The sequence of viral orf10 (UniProt A0A663DJA2, NCBI YP_009725255.1) was downloaded separately from UniProt,
also on June 3, 2020.

