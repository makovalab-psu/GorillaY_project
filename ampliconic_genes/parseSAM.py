import pysam
import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
import sys

def getAssemblyId(assembly):
     if ("frozen" in assembly): 
          return "FRHYB"
     if ("cleaned" in assembly): 
          return "CLHYB"
     if ("assembly15kb" in assembly): 
          return "U15k"
     if ("assembly11kb" in assembly): 
          return "U11k"
     if ("assembly9kb" in assembly): 
          return "U9k"
     if ("assembly7kb" in assembly): 
          return "U7k"
     if ("assembly5kb" in assembly): 
          return "U5k"
     if ("automaticThresholdkb" in assembly): 
          return "UAT"

#python parseSAM.py samfile fastafile
samF=sys.argv[1]
fastaF=sys.argv[2]
geneF=sys.argv[3]
outputfileExtractedGaF=(samF + "_extractedGenesAndFlanks.fasta")
outputfileExtractedG=(samF + "_extractedGenes.fasta")
outputfileReference=(samF + "_reference.txt")
outputfileLinks=(samF + "_links.txt")


samfile = pysam.AlignmentFile(samF, "rb")
fastafile = pysam.Fastafile(fastaF)
gene_dict = SeqIO.to_dict(SeqIO.parse(geneF, "fasta"))

unitigsToOutput=[]

GaF = open(outputfileExtractedGaF, 'w')
G = open(outputfileExtractedG, 'w')
fl = open(outputfileLinks, 'w')

for read in samfile:
     #print read
     if (read.is_unmapped==False):
          genomic_alignments=read.get_blocks()
          #print genomic_alignments
          #rightmost_genomic_coordinate=(genomic_alignments[-1])[1]
          #print rightmost_genomic_coordinate
          #print read.flag
          #print read.cigarstring
          read_id=read.reference_id
          #print read_id
          contig_name=samfile.getrname(read_id)

          array=(read.query_name).split('|')
          orientation_in_human=array[3]
          reversed_in_human=(orientation_in_human=="-1")
          gene_name=array[4].split("/")[0]

          human_start=array[1]
          human_end=array[2]

          parsed_name=(gene_name + " " + human_start + " " + human_end + " " + orientation_in_human)

          start_coord=read.reference_start
          end_coord=read.reference_end
          s=fastafile.fetch(reference=contig_name,start=start_coord,end=end_coord)

          start_coord=max(0,read.reference_start-1000)
          end_coord=read.reference_start
          uf=fastafile.fetch(reference=contig_name,start=start_coord,end=end_coord)

          start_coord=read.reference_end
          end_coord=min(read.reference_end+1000,fastafile.get_reference_length(contig_name)+1)
          df=fastafile.fetch(reference=contig_name,start=start_coord,end=end_coord)
          
          my_seq=SeqRecord(Seq(s,generic_dna))
          upstream_flank=SeqRecord(Seq((uf+"F"),generic_dna))
          downstream_flank=SeqRecord(Seq(("F"+df),generic_dna))

          seq_and_flanks=upstream_flank+my_seq+downstream_flank
          seq_without_flanks=my_seq

          #print (read.query_alignment_length)
          #print (len(gene_dict[read.query_name.split("/")[0]]))
          #CHECK!!!!
          proportion=float(read.query_alignment_length)/(len(gene_dict[read.query_name.split("/")[0]]))
          threshold_lower=0.60
          threshold_upper=1.40

          #print gene_name
          #print read.query_alignment_length
          #print len(gene_dict[read.query_name.split("/")[0]])
          print proportion
          #print read.cigarstring

          assembly_id=getAssemblyId(samF)

          #add assembly name to the header
          contig_id=contig_name.replace("|quiver","")
          contig_id=contig_id.replace("Contig",assembly_id).replace("unitig",assembly_id)

          if (proportion>=threshold_lower and proportion<=threshold_upper):
               #print ("proportion: " + str(proportion))

               #unitig should be included in the reference
               if (contig_name not in unitigsToOutput):
                    #add to the list
                    unitigsToOutput.append(contig_name)

               if (read.is_reverse==False):
                    #reversed in the assembly
                    output=">" + gene_name + "_" + contig_id + "_" + str(read.reference_start+1) + "_" + str(read.reference_end) + "_" + orientation_in_human + "_ASIS\n"
                    GaF.write(output)
                    G.write(output)

                    links = ("gf" + contig_name.replace("|quiver","") + "\t" + str(read.reference_start+1) + "\t" + str(read.reference_end) + "\tgfhumanY\t" + human_start + "\t" + human_end)
                    print (links + "\t" + gene_name)
                    fl.write(links+"\n")
               else:
                    #REVERSE COMPLEMENT
                    seq_and_flanks=(seq_and_flanks.reverse_complement())
                    seq_without_flanks=(seq_without_flanks.reverse_complement())

                    output=">" + gene_name + "_" + contig_id + "_" + str(read.reference_start+1) + "_" + str(read.reference_end) + "_" + orientation_in_human + "_REVC\n"
                    GaF.write(output)
                    G.write(output)

                    links=("gf" + contig_name.replace("|quiver","") + "\t" + str(read.reference_start+1) + "\t" + str(read.reference_end) + "\tgfhumanY\t" + human_start + "\t" + human_end)
                    print (links + "\t" + gene_name)
                    fl.write(links+"\n")
               GaF.write(str(seq_and_flanks.seq)+"\n")
               G.write(str(seq_without_flanks.seq)+"\n")
samfile.close()
GaF.close()
fl.close()

f = open(outputfileReference, 'w')

#output reference file for circos
for i in unitigsToOutput:
     contig=i.replace("|quiver","")
     seq=fastafile.fetch(reference=i)
     print ("chr\t-\tgf"+ contig + "\t" + contig + "\t1\t" + str(len(seq)+1) + "\tchr" + contig)
     #print seq
     f.write("chr\t-\tgf"+ contig + "\t" + contig + "\t1\t" + str(len(seq)+1) + "\tchr" + contig + "\n")

#add human as a reference
f.write("chr\t-\tgfhumanY\thumanY\t1\t57227414\tchrhumanY")
f.close()










