import pysam
import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from collections import defaultdict
import sys

#python spanned_by_long_reads
samF=sys.argv[1]
samfile = pysam.AlignmentFile(samF, "rb")

reference_sequence_confidence=defaultdict(int)
arm1_support=defaultdict(int)
arm2_support=defaultdict(int)

pillow=500 #how many bp needs to read go into arms

for read in samfile:
     #print read
     if (read.is_unmapped==False):
          #get_overlap(
          #return number of aligned bases of read overlapping the interval start and end on the reference sequence.

          read_id=read.reference_id
          reference_name=samfile.getrname(read_id) #full_Contig14:159277-165540_spacer:162256-162617
          array_spacer=reference_name.split(':')[2]
          array_start=reference_name.split(':')[1]

          palindrome_start=int(array_start.split('-')[0])
          spacer_start=int(array_spacer.split('-')[0])-palindrome_start-pillow #pillow makes sure that read captures not only whole spacer, but also part of arm
          spacer_end=int(array_spacer.split('-')[1])-palindrome_start+pillow #pillow makes sure that read captures not only whole spacer, but also part of arm
 

          if (read.reference_start < spacer_start):
               arm1_support[reference_name]=arm1_support[reference_name]+1

          if (read.reference_end > spacer_end):
               arm2_support[reference_name]=arm2_support[reference_name]+1

          if ((read.reference_start < spacer_start) and (read.reference_end > spacer_end)):
               #pacbio support for given palindrome - pacbio reads span spacer
               reference_sequence_confidence[reference_name]=reference_sequence_confidence[reference_name]+1
               print "SUPPORTED"
               print ("READ POSITIONS " + str(read.reference_start) + " " + str(read.reference_end) + " <-> SPACER POSITIONS " + str(spacer_start) + " " + str(spacer_end) + " " + str(reference_name)) 

          #else:
               #print "NOT SUPPORTED"

samfile.close()

#print ("ARM1 support")
#for k,v in arm1_support.items():
#    print k, 'corresponds to', v

#print ("ARM2 support")
#for k,v in arm2_support.items():
#    print k, 'corresponds to', v


print ("ARM1 and ARM2 support: " + str(len(reference_sequence_confidence.items())) + " cases." )
for k,v in reference_sequence_confidence.items():
    print k, 'supported by', v







