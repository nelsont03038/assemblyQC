#!/usr/bin/python

############################################################################
# assemblyQC.py
# written by Thomas Nelson, July 15, 2013 
# thomas dot c dot nelson at gmail dot com
# This script requires the Biopython module and Python version 2.x
# Get the Biopython module here here: http://biopython.org/wiki/Download
# This script will calculate common assembly QC parameters
############################################################################

import sys
from Bio import SeqIO

#usage statement
if len(sys.argv) < 2:
      print "Usage: python assemblyQC.py [ FASTA file name ]"
      sys.exit()

#global variables
num_contigs = 0
contig_len = []

#count the number of contigs and make a sorted list of contig lengths
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    num_contigs += 1
    contig_len.append(len(seq_record.seq))
contig_len.sort()

#get the total size of the "genome"
total_len = sum(contig_len)

#get the mean of the contig lengths
ave_len = total_len / num_contigs

#get the min and max contig lengths
min_len = min(contig_len)
max_len = max(contig_len)

#calculate the median contig length
def median(contig_len):
    if not num_contigs % 2:
        return (contig_len[num_contigs / 2] + contig_len[num_contigs / 2 - 1]) / 2
    return contig_len[num_contigs / 2]
median_len = median(contig_len)

#calculate the NXX statistic (can be any percent, for example the N50)
def nxx(contig_len, fraction):
    cum_len = 0
    for length in contig_len:
        cum_len += length
        if cum_len >= total_len * fraction:
            return length
#change the second argument value to whatever fraction you want to calculate
n50 = nxx(contig_len, 0.5)

#print the output (change the print statements and it should work with Python 3.x)
print "Number of contigs = ", num_contigs
print "Total length = ", total_len
print "Mean contig length = ", ave_len
print "Median contig length = ", median_len
print "Min contig size = ", min_len
print "Max contig size = ", max_len
print "N50 length = ", n50
