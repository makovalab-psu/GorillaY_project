#!/bin/sh

#EDITED


# This is the main RecoverY script that has four parts 
# part1 : run DSK on raw/trimmed reads
# part2	: choose a threshold for the table and build Y k-mer table
# part3 : retrieve Y specific reads of R1
# part4 : retrieve mates for each R1



R1_fsY_reads= $1

R2_fsY_reads= $2



#run DSK on R1 file of raw/trimmed reads, requires dsk2.0.2 or higher 

dsk -file R1_fsY_reads -kmer-size 25 -out R1_fsY_dsk25 -verbose 0

h5dump -y -d dsk/histogram R1_fsY_dsk25.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - >
 R1_fsY.dsk25.histo

dsk2ascii -file R1_fsY_dsk25 -out R1_fsY_dsk25.parse_results  -verbose 0



#pick Threshold
Threshold = 50

awk '$2 > Threshold {print 0}' > post_Threshold_reference_table.fastq

Y_reference_table='./post_Threshold_reference_table.fastq'



echo 'Pointing to all the right places, now classify a subset of R1 reads as Y-chr'

python classify_as_Y_chr.py $R1_fsY_reads $Y_reference_table post_RecoverY/RecY_pairs.A.fastq

echo 'Done with classify, now retrieve the mate for each read'

python find_mates.py post_RecoverY/RecY_pairs.A.fastq $R2_fsY_reads post_RecoverY/RecY_pairs.B.fastq

echo 'Done with RecoverY'

