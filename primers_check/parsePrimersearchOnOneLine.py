import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from difflib import SequenceMatcher
#python parsePrimersearchOnOneLine.py primersearch assembly_fasta

def similar(a, b):
	return SequenceMatcher(None, a, b).ratio()

primersearchF=sys.argv[1]
assemblyF=sys.argv[2]
primersF=sys.argv[3]

af=SeqIO.to_dict(SeqIO.parse(assemblyF, "fasta"))
outfile=primersearchF.replace("_temp","")+"_extracted"

#read primers to the dictionary
F_primer_dict=dict()
prf = open(primersF, 'r')
for line in prf:
	line=line.rstrip()
	array=line.split('\t')
	F_primer_dict[array[0]]=array[1]

prf.close()

pfw = open(outfile, 'w')
pf = open(primersearchF, 'r')
for line in pf:
	line=line.rstrip()
	gene_name=line.split(' ')[2]

	amplimers=line.split("Amplimer")

	for i in range(1,len(amplimers)-1,2):

		amplimer=amplimers[i]
		amplimer_length=int(amplimers[i+1].split(" ")[2])

		array=amplimer.split(" ")
		forward_primer_sequence=F_primer_dict[gene_name]

		contig=array[3]
		amplicon_sequence=af[contig].seq
		contig_length=len(af[contig].seq)

		start=int(array[9])
		rev_primer_hit=int(array[18].replace("[","").replace("]",""))
		#print rev_primer_hit
		end=contig_length-rev_primer_hit+1

		min_start=min(start,end)
		max_end=max(start,end)

		print ("gf" + contig + "\t" + str(min_start) + "\t" + str(max_end) + "\t" + gene_name)

		orientation="ASIS"

		#upstream=amplicon_sequence[(min_start-1002):(min_start-2)]
		amplicon=amplicon_sequence[(min_start-1):max_end] #coordinates match samtools faidx
		#downstream=amplicon_sequence[(max_end+1):(max_end+1001)]

		#sequence=upstream+"E"+amplicon+"E"+downstream
		sequence=amplicon
		#print sequence

		#should be reverse complemented?
		in_contig=sequence.lower()
		in_primer=forward_primer_sequence.lower()

		#print in_contig
		#print in_primer
		similarity=similar(in_contig[:len(in_primer)],in_primer)
		#print similarity

		if (similarity<0.9):
			sequence=sequence.reverse_complement()
			orientation="REVC"

		#print sequence
		#print "\n"

		if (len(sequence)<=1000):
			pfw.write(">" + gene_name + "_" + contig + "_" + str(min_start) + "_" + str(max_end) + "_" + orientation + "\n") #NOTE: we don't know strand
			pfw.write(str(sequence) + "\n")

pfw.close()
pf.close()