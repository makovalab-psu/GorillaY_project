#!/bin/bash -
#folder with sequences, folder with assemblies
mkdir alignments figures;

#for all sequences and all assemblies
for seq in ${1}/*.fasta; do 

	#define name variables 
	seq_name=`echo $(basename $seq) | sed 's/.fasta//g'`;
	rm -f ${seq_name}_merged.txt; #remove merged file if existed

	for as in ${2}/*.fasta; do 

		assembly_name=`echo $(basename $as) | sed 's/.fasta//g'`;
		outputfile=${seq_name}_${assembly_name}; 
		sam=alignments/${outputfile}.sam;

		echo "seq:" $seq $seq_name;
		echo "as:" $as $assembly_name;
		echo "out:" $outputfile;
		echo "sam:" $sam;

		bwa mem -k 5 $as $seq >$sam;
		python calculateAlignment.py $seq $sam;

		echo -n "gene_name length "$assembly_name" " >>${seq_name}_merged.txt;
	done; 

	#combine results from all assemblies into single file
	echo "" >>${seq_name}_merged.txt; #add newline after header
	paste -d ' ' alignments/*$seq_name*toPlot.txt  >>${seq_name}_merged.txt;
done;

#plot in R
for file in *merged.txt; do
	echo "plotting:" $file; perl -pi -w -e 's/ $//g;' $file;
	Rscript gene_screening_heatmapByGenes.R $file;
	#Rscript gene_screening_heatmapByAssembly.R $file;
done;