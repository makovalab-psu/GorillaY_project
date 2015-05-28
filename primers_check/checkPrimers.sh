#!/bin/bash
#checkPrimers.sh primers_file mismatch_percentage
set -e
set -x
primers_file=$1;
mismatch=$2

function search_primers {
	primers_name=$(basename $primers_file)
	fasta_name=$(basename $1)
	primersearch -mismatchpercent $mismatch -infile $primers_file -seqall $1 -outfile primersearch.${primers_name}.${fasta_name}.${mismatch} &
} 

#check for presence of primers in female, humanY and chimpY and assembly itself

#female
search_primers /galaxy/home/biomonika/data/ref_gorilla/Gorilla_gorilla.gorGor3.1.74.dna.toplevel.fa;

#Y chromosomes
search_primers /galaxy/home/biomonika/data/ref_Ys/chrY_human.fa;
search_primers /galaxy/home/biomonika/data/ref_Ys/chrY_chimp.fa;
search_primers /galaxy/home/biomonika/data/ref_Ys/chrY_rhesus.fa;

#assemblies
search_primers /nfs/brubeck.bx.psu.edu/scratch4/monika/summary_of_assemblies/Filtered/4.pbjelly_filt.fasta_1000
wait;