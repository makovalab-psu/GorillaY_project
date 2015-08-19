#!/bin/sh

#requires dsk2.0.2 or higher 
#run as ./run_dsk.sh FASTA_file_to_be_kmerized

if [ $# -ne 1 ]
then
	echo "Usage: $0 <FASTA_file_to_be_kmerized>"
	exit 1
fi


R1_fsY_reads=$1

echo $R1_fsY_reads

dsk -file $R1_fsY_reads -abundance-min 1 -kmer-size 11 -out R1_fsY_dsk -verbose 0

h5dump -y -d dsk/histogram R1_fsY_dsk.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > R1_fsY.dsk.histo

dsk2ascii -file R1_fsY_dsk -out R1_fsY_dsk.parse_results  -verbose 0