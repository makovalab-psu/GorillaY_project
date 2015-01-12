#!/usr/local/bin/python
# fill gaps in scaffolds based on reference (e.g. human) 
#python walkAssembly.py genesHumanBiomartY.list genes_on_assembly.sam genesHumanBiomartY.fasta assembly_in_fasta
import sys
import test
import os
import re
import subprocess
from Bio import SeqIO

class GeneMappings:
    """Class that for every gene stores where it maps"""
    def __init__(self, name, start, end, orientation, mappings, number_of_mappings, flagstat, mapping_positions):
        self.name = name
        self.start = start
        self.end = end
        self.orientation = orientation
        self.mappings = mappings
        self.number_of_mappings=number_of_mappings
        self.flagstat = flagstat
        self.mapping_positions = mapping_positions

def getMappingsForGene(gene): 
    command = "cat " + sam_file + " | grep -v \"@\" | grep \"" + gene + "\" | cut -f3 "
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True)
    #Launch shell command:
    gene_mappings = process.communicate()[0]

    command2 = "cat " + sam_file + " | grep -v \"@\" | grep \"" + gene + "\" | cut -f2 | awk '{s+=$1} END {print s}'"
    process2 = subprocess.Popen(command2, stdout=subprocess.PIPE, stderr=None, shell=True)
    #Launch shell command:    
    cumulative_flag_score = process2.communicate()[0]

    command3 = "cat " + sam_file + " | grep -v \"@\" | grep \"" + gene + "\" | cut -f4 "
    process3 = subprocess.Popen(command3, stdout=subprocess.PIPE, stderr=None, shell=True)
    #Launch shell command:    
    mapping_positions = process3.communicate()[0]

    return (gene_mappings + "#" + cumulative_flag_score + "#" + mapping_positions)
    #return command

gene_file = sys.argv[1]
sam_file = sys.argv[2]
fasta_file = sys.argv[3]
assembly_file = sys.argv[4]

print ("gene file: " + gene_file)
print ("sam file: " + sam_file)

output_file = sam_file + "_gapFilling"

toWalk=open(output_file+".txt", "w")
toWalkFasta=open(output_file+".fasta", "w")
toWalkStat=open(output_file+".stat.txt", "w")

handle = open(fasta_file, "rU")
record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()

handle = open(assembly_file, "rU")
record_assembly = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()


#for every gene in order of humanY
with open(gene_file, "r") as outer:
        ROUND=0
        for outerline in outer:
            ROUND=ROUND+1
            print ("ROUND: " + str(ROUND))

            array=outerline.split('\t')
            gene=array[4]
            gene=gene.rstrip('\n')

            gene_array=getMappingsForGene(gene)
            gene_array=gene_array.replace('\n', ';')
            #divide into mappings and flagstat
            gene_array=gene_array.split('#')

            geneOuter=GeneMappings(name=gene, start=int(array[1]), end=int(array[2]), orientation=int(array[3]), mappings=gene_array[0], number_of_mappings=gene_array[0].count(';'), flagstat=int(gene_array[1].replace(';', '')), mapping_positions=gene_array[2].replace(';', ''))


            with open(gene_file, "r") as inner:
                for innerline in inner:
                    array=innerline.split('\t')
                    gene=array[4]
                    gene=gene.rstrip('\n')

                    gene_array=getMappingsForGene(gene)
                    gene_array=gene_array.replace('\n', ';')
                    #divide into mappings and flagstat
                    gene_array=gene_array.split('#')
                    geneInner=GeneMappings(name=gene, start=int(array[1]), end=int(array[2]), orientation=int(array[3]), mappings=gene_array[0], number_of_mappings=gene_array[0].count(';'), flagstat=int(gene_array[1].replace(';', '')), mapping_positions=gene_array[2].replace(';', ''))

                    #look up contigs in fasta 
                    outer_contig=geneOuter.mappings.replace(';', '');
                    inner_contig=geneInner.mappings.replace(';', '');

                    if (geneOuter.number_of_mappings + geneInner.number_of_mappings == 2 ) and (geneOuter.mappings != geneInner.mappings) and (outer_contig != "*") and (inner_contig != "*") and (geneOuter.end <= geneInner.start):

                        print ("outer gene: " + geneOuter.name + " innner gene: " + geneInner.name)

                        #retrieve assembly sequences
                        outer_record=record_assembly[outer_contig] #use any record ID
                        inner_record=record_assembly[inner_contig]

                        #calculate distance between contigs
                        #first calculate distance between genes

                        #lengths of contigs
                        length_of_outer_contig=int(len(outer_record.seq))
                        length_of_inner_contig=int(len(inner_record.seq))

                        PCRdistanceBetweenGenes=geneInner.start - geneOuter.end #this distance is over-estimation of real distance, we need to subtract distances of the genes to the end of contigs
                        print ("gene distance estimate: " + str(PCRdistanceBetweenGenes))

                        if (PCRdistanceBetweenGenes>500000):
                            print ("current gene too distant from outer gene, trying new outer gene; break;")
                            break;

                        #subtract distance of geneA to the end of contig
                        distance_to_the_end_of_first_contig=(length_of_outer_contig-int(geneOuter.mapping_positions))
                        PCRdistanceEstimate=PCRdistanceBetweenGenes-distance_to_the_end_of_first_contig

                        #substract distance of geneB to the beginning of contig
                        distance_to_the_start_of_gene_at_second_contig=int(geneInner.mapping_positions)
                        PCRdistanceEstimate=PCRdistanceEstimate-distance_to_the_start_of_gene_at_second_contig

                        print ("distance_to_the_end_of_first_contig: " + str(distance_to_the_end_of_first_contig) + " distance_to_the_start_of_gene_at_second_contig: " + str(distance_to_the_start_of_gene_at_second_contig))

                        print ("PCR distance estimate: " + str(PCRdistanceEstimate))

                        if ((PCRdistanceEstimate<=20000) and (PCRdistanceEstimate>-10000)):
                            print (str(PCRdistanceEstimate) + " PCR product size OK")
                            print outerline
                            print vars(geneOuter)
                            print innerline
                            print vars(geneInner)

                            #check for orientation
                            print (str(geneOuter.orientation) + " " + str(geneInner.orientation))
                            orientation_score=geneInner.orientation+geneOuter.orientation
                            print ("orientation_score: " + str(orientation_score))

                            print ("outer flagstat: " + str(geneOuter.flagstat) + " flagstat: " + str(geneInner.flagstat))

                            #reverse complement if needed
                            if ((geneOuter.orientation==1) and (geneInner.orientation==1)):
                                human_orientation="sense sense"
                                print human_orientation
                                #both in sense orientation in human, so both need to be in sense orientation in contigs
                                if (geneOuter.flagstat==16):
                                    outer_record.seq=(outer_record.seq).reverse_complement()
                                    print "reverse complementing outer gene"
                                if (geneInner.flagstat==16):
                                    record.seq=(inner_record.seq).reverse_complement()
                                    print "reverse complementing gene"

                            if ((geneOuter.orientation==-1) and (geneInner.orientation==-1)):
                                human_orientation="antisense sense"
                                print human_orientation
                                #both in reverse orientation in human, so both need to be in reverse orientation in contigs
                                if (geneOuter.flagstat==0):
                                    outer_record.seq=(outer_record.seq).reverse_complement()
                                    print "reverse complementing outer gene"
                                if (geneInner.flagstat==0):
                                    record.seq=(inner_record.seq).reverse_complement()
                                    print "reverse complementing gene"

                            if ((geneOuter.orientation==1) and (geneInner.orientation==-1)):
                                human_orientation="sense antisense"
                                print human_orientation
                                if (geneOuter.flagstat==16):
                                    outer_record.seq=(outer_record.seq).reverse_complement()
                                    print "reverse complementing outer gene"
                                if (geneInner.flagstat==0):
                                    record.seq=(inner_record.seq).reverse_complement()
                                    print "reverse complementing gene"

                            if ((geneOuter.orientation==-1) and (geneInner.orientation==1)):
                                human_orientation="antisense sense"
                                print human_orientation
                                #both in sense orientation in human, so both need to be in sense orientation in contigs
                                if (geneOuter.flagstat==0):
                                    outer_record.seq=(outer_record.seq).reverse_complement()
                                if (geneInner.flagstat==16):
                                    record.seq=(inner_record.seq).reverse_complement()

                            #retrieve updated sequences
                            outer_fasta=("\n>" + outer_record.id + "\n" + outer_record.seq)
                            inner_fasta=("\n>" + inner_record.id + "\n" + inner_record.seq)

                            candidate_report="\nCANDIDATE " + geneOuter.name + " maps to " + geneOuter.mappings + " " + geneInner.name + " maps to " + geneInner.mappings + " distance_between_genes is: " + str(PCRdistanceBetweenGenes) + "; PCR product size is: " + str(PCRdistanceEstimate) + "; distance_from_gene_to_end_of_first_contig: " + str(distance_to_the_end_of_first_contig) + " distance_to_start_of_gene_at_second_contig: " + str(distance_to_the_start_of_gene_at_second_contig) + "; length_of_first_contig: " + str(length_of_outer_contig) + "; length_of_second_contig: " + str(length_of_inner_contig) + "; human_orientation: "  + human_orientation;
                            
                            print candidate_report; 
                            toWalk.write(str(candidate_report) + str(outer_fasta) + str(inner_fasta))

                            #retrieve genes
                            outer_record_genes=record_dict[geneOuter.name] #use any record ID
                            inner_record_genes=record_dict[geneInner.name]
                        
                            outer_fasta_gene=(">" + outer_record_genes.id + "\n" + outer_record_genes.seq)
                            inner_fasta_gene=(">" + inner_record_genes.id + "\n" + inner_record_genes.seq)
                            
                            gene_report=("\n" + outer_fasta_gene + "\n" + inner_fasta_gene + "\n----------")
                            print (gene_report)
                            toWalk.write(str(gene_report))
                            toWalkFasta.write(str(outer_fasta) + str(inner_fasta))
                            toWalkStat.write(str(candidate_report) + "\n")
                        else:
                            print (str(PCRdistanceEstimate) + " PCR product size NOK")

                    #else:
                        #print (" NOT CANDIDATE")

# Close opened file
toWalk.close()
toWalkFasta.close()
toWalkStat.close()




