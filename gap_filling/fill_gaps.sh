#!/bin/bash
#./fill_gaps.sh assembly_in_fasta genesHumanBiomartY.list genesHumanBiomartY.fasta

{

date >&2
set -x;

assembly_in_fasta=$1
genes_list=$2
genes_fasta=$3
sam=${genes_fasta}_Genes_Mapped_To_${assembly_in_fasta}.sam;
logfile=${genes_fasta}_Genes_Mapped_To_${assembly_in_fasta}.logfile.txt

#map genes on assembly and create sam file

if [ -s $sam ];
then
   echo "File $sam exists. No mapping will be performed."
else
   echo "File $sam does not exist or is empty Mapping will be performed."
   bwa index $assembly_in_fasta;
   bwa mem -t 60 $assembly_in_fasta $genes_fasta >$sam;
fi

wait;
#fill gaps and write result to the file
#python walkAssembly.py genesHumanBiomartY.list genes_on_assembly.sam genesHumanBiomartY.fasta assembly_in_fasta
python PCRdistance.py $genes_list $sam $genes_fasta $assembly_in_fasta;

date >&2

} 2>fill_gaps.log;
