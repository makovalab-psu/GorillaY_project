#!/usr/local/bin/python
# script reads fasta file and sam file and output file with 3 columns - gene name, length, alignment length.
import sys
import os
from Bio import SeqIO
import subprocess

def getAlignmentLengthForGene(gene): 
    command = "cat " + sam_file + "| grep -v \"@\" | grep \"\\b" + gene + "\\b\" | perl -lane '$l = 0; $F[5] =~ s/(\d+)[MX=N]/$l+=$1/eg; print $l' | awk '{ sum+=$1} END {print sum}'";
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True)
    #Launch the shell command:
    output = process.communicate()[0];
    return output;

fasta_file = sys.argv[1]
sam_file = sys.argv[2]
output_file = sam_file + "toPlot.txt";

toPlot=open(output_file, "wb");

#for every sequence in fasta file
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    print(seq_record.id)
    print(len(seq_record))
    #print name of sequence
    toPlot.write(seq_record.id + " ");
    #print length of sequence
    toPlot.write(str(len(seq_record)) + " ");
    #print alignment length
    toPlot.write(str(getAlignmentLengthForGene(seq_record.id)));

# Close opened file
toPlot.close()

