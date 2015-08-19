#!/bin/bash

# ./recovery_main.sh READS_R1.FASTQ READS_R2.FASTQ PRE_THRESHOLD_REF_TBL THRESHOLD_VALUE

# This is the main RecoverY script that has 3 parts 
# part1	: choose a threshold for the table and build Y k-mer table
# part2 : retrieve Y specific reads of R1
# part3 : retrieve mates for each R1

if [ $# -ne 4 ]
then
	echo "Usage: $0 READS_R1.FASTQ READS_R2.FASTQ PRE_THRESHOLD_REF_TBL THRESHOLD_VALUE"
	exit 1
fi


R1_fsY_reads=$1

R2_fsY_reads=$2

pre_Threshold_reference_table=$3

Threshold=$4

mkdir post_RecoverY

# -------> START part1 <------

awk -v t=$Threshold '$2>t {print $0}' $pre_Threshold_reference_table > tables/post_Threshold_reference_table

Y_reference_table='tables/post_Threshold_reference_table'


# -------> END OF part1 <------



# -------> START part2 <------

echo 'Pointing to all the right places, now classify a subset of R1 reads as Y-chr'

python py_scripts/classify_as_Y_chr.py $R1_fsY_reads $Y_reference_table post_RecoverY/RecY_pairs.A.fastq

echo 'Done with classify, now retrieve the mate for each read'

# -------> END OF part2 <------



# -------> START part3 <------

python py_scripts/find_mates.py post_RecoverY/RecY_pairs.A.fastq $R2_fsY_reads post_RecoverY/RecY_pairs.B.fastq

echo 'Done with RecoverY'
