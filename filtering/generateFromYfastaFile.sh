#!/bin/bash
#./generateFromYfastaFile.sh folder


{
date >&2
set -x;
shopt -s extglob; 
folder=$1;
path_to_reference="/scratch/assemblies/Gorilla+YY";
hash blasr 2>/dev/null || { echo "I require blasr in path but it's not installed.  Aborting."; exit 1; }
hash seqtk 2>/dev/null || { echo "I require seqtk in path but it's not installed.  Aborting."; exit 1; }

#for all fasta files in folder
for file in *.f*a; do 
	echo "File being processed: " $file; 

	if [[ $file =~ .*fromY.* ]]
	then
		echo "FromY fasta file won't be processed.";
	else
		#check if file has already been mapped to Gorilla+YY (female gorilla AND human Y AND chimp Y) with blasr and bestn=1, only best hit gets reported
		if [ -f ${file}_mapped_Gorila+YY ];
		then
			echo "File ${file}_mapped_Gorila+YY exists.";
		else
			echo "File ${file}_mapped_Gorila+YY does not exist. Blasr will be used to map contigs to Gorilla+YY";
			blasr -bestn 1 -nproc 16 -minPctIdentity 70 -m 4 -sa ${path_to_reference}.sa $file $path_to_reference >${file}_mapped_Gorila+YY;
		fi

		#create contamination file - contigs with best autosomal hit
		cat ${file}_mapped_Gorila+YY | grep -v "@" | awk '{if ($2!="X" && $2!="chimpY" && $2!="humanY" && $2!="unplaced") {print $0} }' |  cut -d' ' -f1 | cut -d '/' -f1 | sort >contamination; 
	
		#create list of all contigs
		cat ${file} | grep ">" | sed 's/>//' | cut -d' ' -f1 | sort >all; 
	
		#filter out contamination, keep Y-homologs and newly discovered sequences
		comm -13 contamination all >trueY; 
	
		#generate fromY fasta file
		seqtk subseq $file trueY >${file}_fromY.fasta; 
	
		#count how many contigs are in original and fromY files
		cat $file | grep ">" | wc -l; 
		cat ${file}_fromY.fasta | grep ">" | wc -l; 

		#remove files
		#rm contamination all trueY;
	fi
done; 
date >&2

} 2>generateFromYfastaFile.log;
