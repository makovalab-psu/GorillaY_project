import sys
import test
import os
import re
import subprocess
import os.path
from Bio import SeqIO
#extract_palindromes.py assembly_file

assembly_file=sys.argv[1]
assembly_name=assembly_file.replace(".fasta","")
palindrome_output=(assembly_file + "_assemblyNewPalindromes.fasta")
statistics_output=(assembly_file + "_assemblyNewPalindromes_statistics.txt")
RMfile=palindrome_output+".masked"
RMparsed="RM_parsed.txt"
bestHitfile="blast.xml"
bestHitfileparsed="blast.xml_parsed.txt"

handle = open(assembly_file, "rU")
record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()

fasta_output = open(palindrome_output,"w")
stat_output = open(statistics_output,"w")

#BLASR output
def writeToBlasrFile():
    command = "blasr -bestn 2 -nproc 60 -m 4 -out BLASRfemale_aligned.txt -unaligned BLASRfemale_unaligned.txt " + palindrome_output + " /galaxy/home/biomonika/data/ref_gorilla/Gorilla_gorilla.gorGor3.1.74.dna.toplevel.fa -sa /galaxy/home/biomonika/data/ref_gorilla/gorGor3.sa"
    #command = "blasr -bestn 2 -nproc 60 -m 4 " + palindrome_output + " /galaxy/home/biomonika/data/ref_gorilla/Gorilla_gorilla.gorGor3.1.74.dna.toplevel.fa >" + BLASRfemale_aligned.txt
    print command
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True) 
    process.communicate()[0] #Launch shell command

def prepareBlasrOutput():
    if (os.path.isfile('BLASRfemale_aligned.txt') and os.path.getsize('BLASRfemale_aligned.txt') > 0): #blast file exists and is not empty
        with open('BLASRfemale_aligned.txt', 'r') as content_file:
            BlasrFile = content_file.read()
            if 'ERROR' in BlasrFile:
                content_file.close()
                writeToBlasrFile()
            else:
                print ("blasr mapping to female exists, blasr won't be run: " + "BLASRfemale_aligned.txt")
    else: 
        print ("blasr mapping DOES NOT exist. Blasr will be run: " + "BLASRfemale_aligned.txt")
        writeToBlasrFile()


#BLAST output
#blastn seq versus seq.r to identify inverted repeat, hence potential palindromes
blast_file=(assembly_name+"__palindromes_defWordSize.txt")
if (os.path.isfile(blast_file) and os.path.getsize(blast_file) > 0): #blast file exists and is not empty
    print ("blastfile exists, blastn won't be run:" + blast_file)
else:
    if (os.path.isfile(assembly_name+".nsq")):
        print "Blast database already exists."
    else:
        print "Creating blast database."
        command = "makeblastdb -in " + assembly_file + " -out " + assembly_name + " -dbtype nucl"
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True) 
        process.communicate()[0] #Launch shell command

    #blast assembly against assembly reverse complement
    query=(assembly_file+".r") #reverse complemented assembly is query
    if (os.path.isfile(query) and os.path.getsize(query) > 0): #reversed complemented file exists and is not empty
        print ("Reverse complemented assembly exists, won't be created:" + query)
    else:
        print ("Creating reverse complemented assembly sequence: " + query)
        command="seqtk seq -r " + assembly_file +" >" + query
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True) 
        process.communicate()[0] #Launch shell command

    command="blastn -num_threads 30 -db " + assembly_name + " -query " + query + " -outfmt 7 -dust no -perc_identity 95 -strand plus | sort -gk1 >" + blast_file
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True) 
    process.communicate()[0] #Launch shell command

blast_file_parsed=(blast_file + "_parsed.txt")

#sort blast file to keep only reciprocal blast hits for the same contigs/scaffolds
command = "cat " + blast_file + " | grep -v \"#\" | awk '{if ($1==$2) print;}' | sort -gk1,1 -k4,4 >" + blast_file_parsed
process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True) 
process.communicate()[0] #Launch shell command

def getBestBlastHit(scaf_name):
    if (os.path.isfile(bestHitfileparsed) and os.path.getsize(bestHitfileparsed) > 0): #RM file exists and is not empty
        with open(bestHitfileparsed,"r") as BH_file:
            for line in BH_file:
                if scaf_name in line:
                    forbidden = ['F1','A1','S1','A2','F2']
                    try:
                        nextline=next(BH_file)
                        if any(forb in nextline for forb in forbidden):  #check if blast hit exists for this entry
                            return "none"
                        else:
                            return nextline.rstrip()
                    except StopIteration:
                        break

def getRMscore(scaf_name):
    if (os.path.isfile(RMparsed) and os.path.getsize(RMparsed) > 0): #RM file exists and is not empty
        with open(RMparsed, 'r+') as RM_file:
            RM = RM_file.read()
            for line in RM.splitlines():
                #print line
                if scaf_name in line:
                    array=line.split(' ')
                    return "%RM: " + array[1]

def getHitAndLOA(scaf_name):
    hit=0
    ret1=[0,0,0,0]
    ret2=[0,0,0,0]

    with open('BLASRfemale_aligned.txt', 'a+') as content_file:
        BlasrFile = content_file.read()

        for line in BlasrFile.splitlines():
            if scaf_name in line:
                array=line.split(' ')
                l=(float(array[6])-float(array[5]))/(float(array[7])) #length of alignment
                length=int(round(l*100))

                if (hit==0):
                    #first hit found
                    ret1=[array[1],length,array[9],array[10]]  
                    #print (scaf_name + " FOUND IN " + line)
                    #print ("HIT: " + array[1] + " L: " + str(length) + " FEMALE_POS: " + array[9] + " " + array[10] )

                if (hit==1):
                    #second hit found
                    ret2=[array[1],length,array[9],array[10]] #chrname,%length_of_alignment,female_start,female_end
                    #print (scaf_name + " FOUND IN " + line)
                    #print ("HIT: " + array[1] + " L: " + str(length) + " FEMALE_POS: " + array[9] + " " + array[10] )
                hit+=1
        return [ret1,ret2]

def getDistInFemale(first, second):
    hitIndexi=[0,1] #we are interested in two best hits to female genome
    hitIndexj=[0,1]

    first_name=first.partition("_")[0]
    second_name=second.partition("_")[0]

    #there are 4 combinations if two best hits for both are considered
    for i in hitIndexi:
        for j in hitIndexj:
            print ('Combining hits ' + str(i+1) + ' and ' + str(j+1) + ':'),
            if ((getHitAndLOA(first)[i][1]>0) and (getHitAndLOA(second)[j][1]>0)): #both alignments are greater than 0 => alignments exist and distance can be calculated
                if ((getHitAndLOA(first)[i][0])==(getHitAndLOA(second)[j][0])):
                    dist = int(getHitAndLOA(second)[i][3]) - int(getHitAndLOA(first)[j][3])
                    print ("distance from " + first_name + " to " + second_name + ": " + str(dist))
                else: 
                    print ("distance from " + first_name + " to " + second_name + ": not applicable, different chromosomes")
            else:
                print ("distance from " + first_name + " to " + second_name + ": nan, one of the sequences doesn't exists")

def getOverlap(a, b):
    return max(0, min(int(a[1]), int(b[1])) - max(int(a[0]), int(b[0])))

def checkIfPalindrome(previous_line,line,scaf_len): 

    #check if alignment length matches
    alignment_length=(previous_line[3]==line[3])

    #check if not overlapping on plus strand

    overlapping=(getOverlap([previous_line[6],previous_line[7]],[line[6],line[7]])>0)

    #check if middle part makes sense

    #check if left pillow size makes sense
    A=min([int(previous_line[6]),int(previous_line[7]),int(line[6]),int(line[7])])
    B=scaf_len-max([int(previous_line[8]),int(previous_line[9]),int(line[8]),int(line[9])])
    diff=abs(A-B)
    print ("leftPillow; A: " + str(A) + " B: " + str(B) + " diff: " + str(diff))
    leftPillow=(diff<100)

    #check if right pillow size makes sense
    A=scaf_len-max([int(previous_line[6]),int(previous_line[7]),int(line[6]),int(line[7])])
    B=min([int(previous_line[8]),int(previous_line[9]),int(line[8]),int(line[9])])
    diff=abs(A-B)
    print ("rightPillow; A: " + str(A) + " B: " + str(B) + " diff: " + str(diff))
    rightPillow=(diff<100)

    print("scaf_len: " + str(scaf_len))
    print ("alignment length: " + str(alignment_length) + " overlapping: " + str(overlapping) + " leftPillow: " + str(leftPillow) + " rightPillow: " + str(rightPillow))
    conlusion=(alignment_length and (not overlapping) and leftPillow and rightPillow)
    print ("IS_PALINDROME: " + str(conlusion))
    return conlusion

def getInfoForBlastPair(previous_line,line,identity): 

    strand="NAN"
    if ((int(line[8])-int(previous_line[9])+1)>=0): #spacer_length=start-prev_end+1;
        strand="POSITIVE_STRAND"
        prev_start=min(int(previous_line[8]),int(previous_line[9]))
        prev_end=max(int(previous_line[8]),int(previous_line[9]))

        start=min(int(line[8]),int(line[9]))
        end=max(int(line[8]),int(line[9]))
    else:
        strand="NEGATIVE_STRAND"
        start=min(int(previous_line[8]),int(previous_line[9]))
        end=max(int(previous_line[8]),int(previous_line[9]))

        prev_start=min(int(line[8]),int(line[9]))
        prev_end=max(int(line[8]),int(line[9]))


    scaf_name=line[0]
    
    spacer_length=start-prev_end+1;
    
    #print (str(prev_start) + " " + str(prev_end) + " " + str(start) + " " + str(end))
    #print (scaf_name + ":" + str(prev_end) + "-" + str(start) + " spacer_length:" + str(spacer_length))
    
    #print ("previous_line:" + ', '.join(previous_line))
    #print ("line:" + ', '.join(line))


    #print "PALINDROME:"
    first_flanking_end=prev_start-1
    first_flanking_start=first_flanking_end-299
    first_flanking_start=max(first_flanking_start,1) #maybe it's end of scaffold and therefore no flanking sequence can be extracted
    first_flanking_end=max(first_flanking_end,1) #maybe it's end of scaffold and therefore no flanking sequence can be extracted

    second_flanking_start=end+1
    second_flanking_end=second_flanking_start+299 #if out of the range, won't be extracted with faidx [todo:fix, read sequences to dictionary and get length]

    F1=("F1_" + scaf_name + ":" + str(min(first_flanking_start,first_flanking_end)) + "-" + str(max(first_flanking_start,first_flanking_end)))
    A1=("A1_" + scaf_name + ":" + str(min(prev_start,prev_end)) + "-" + str(max(prev_start,prev_end)))
    S1=("S1_" + scaf_name + ":" + str(min(prev_end+1,start-1)) + "-" + str(max(prev_end+1,start-1)))
    A2=("A2_" + scaf_name + ":" + str(min(start,end)) + "-" + str(max(start,end)))
    F2=("F2_" + scaf_name + ":" + str(min(second_flanking_start,second_flanking_end)) + "-" + str(max(second_flanking_start,second_flanking_end)))

    print "#chrname,%length_of_alignment,female_start,female_end"
    print (F1 + " length:" + str(first_flanking_end-first_flanking_start+1) + "\t" + strand + "\t" + str(getRMscore(F1)) + "\t" + str(getHitAndLOA(F1)) + " BestBlastHit: " + str(getBestBlastHit(F1)))
    print (A1 + " length:" + str(prev_end-prev_start+1) + "\t" + strand + "\t" + str(getRMscore(A1)) + "\t" +  str(getHitAndLOA(A1)) + " BestBlastHit: " + str(getBestBlastHit(A1)))
    print (S1 + " length:" + str((start-1)-(prev_end+1)+1) + "\t" + strand + "\t" + str(getRMscore(S1)) + "\t" +  str(getHitAndLOA(S1)) + " BestBlastHit: " + str(getBestBlastHit(S1)))
    print (A2 + " length:" + str(end-start+1)+ "\t" + strand + "\t" + str(getRMscore(A2)) + "\t" +  str(getHitAndLOA(A2)) + " BestBlastHit: " + str(getBestBlastHit(A2)))
    print (F2 + " length:" + str(second_flanking_end-second_flanking_start+1) + "\t" + strand + "\t" + str(getRMscore(F2)) + "\t" +  str(getHitAndLOA(F2)) + " BestBlastHit: " + str(getBestBlastHit(F2)))

    #write down length of arm1
    stat_output.write(A1 + " length: " + str(prev_end-prev_start+1) + " identity: " + identity + "\n")

    #check adjacency
    #A1 to A2
    getDistInFemale(A1,A2)
    print "==="

    #F1 to A1,A2
    getDistInFemale(F1,A1)
    getDistInFemale(F1,A2)
    print "==="

    #S1 to A1,A2
    getDistInFemale(A1,S1)
    getDistInFemale(S1,A2)
    print "==="

    #F2 to A1,A2
    getDistInFemale(A1,F2)
    getDistInFemale(A2,F2)
    print "==="

    #print sequences to the file
    seqEntry=str(record_dict[scaf_name].seq)

    fasta_output.write(">F1_" + F1 + "\n" + seqEntry[(first_flanking_start-1):(first_flanking_end)-1] + "\n") # items start through end-1, for example 0:5 means 5 bp from 1st bp
    fasta_output.write(">A1_" + A1 + "\n" + seqEntry[(prev_start-1):(prev_end-1)] + "\n")
    fasta_output.write(">S1_" + S1 + "\n" + seqEntry[(prev_end):(start-2)] + "\n")
    fasta_output.write(">A2_" + A2 + "\n" + seqEntry[(start-1):(end-1)] + "\n")
    fasta_output.write(">F2_" + F2 + "\n" + seqEntry[(second_flanking_start-1):(second_flanking_end-1)] + "\n")

previous_scaf1=""
previous_line=""

with open(blast_file_parsed, "r") as f:
    palindromes_found=0
    for line in f:
        array=line.split('\t')
        scaf1=array[0]
        scaf_len=len(str(record_dict[scaf1].seq))
        identity=array[2]
        line=array
        #print scaf1
        if(previous_scaf1==scaf1):
            print (previous_line)
            print (line)
            print ("match: " + previous_scaf1 + " AND " + scaf1)
            print ("inverted repeat on same scaffold detected, now checking if palindrome")

            if (checkIfPalindrome(previous_line,line,scaf_len)):
                print ("candidate palindrome")
                palindromes_found=palindromes_found+1
                getInfoForBlastPair(previous_line,line,identity)
            
            print ("\n\n")

        previous_scaf1=scaf1
        previous_line=line

    print ("\nPALINDROMES FOUND: " + str(palindromes_found) + "\n")

fasta_output.close()
stat_output.close()

#run Blasr if not ready by now
prepareBlasrOutput()

print "NOW RUN RepeatMasker -species Primates -x " + palindrome_output + " and copy " + palindrome_output + ".masked to the actual folder.";
print "BLAST " + palindrome_output + " with nucleotide blast and download xml file. Name it blast.xml and move to the actual folder."
raw_input("Press Enter to continue computation...")

if (os.path.isfile(RMfile) and os.path.getsize(RMfile) > 0): #check if RepeatMasker file exists and is not empty
    print "RepeatMasker file exists, checked."
else:
    print "RepeatMasker file doesn't exist. Script will be terminated."
    sys.exit()

if (os.path.isfile(bestHitfile) and os.path.getsize(bestHitfile) > 0): #check if RepeatMasker file exists and is not empty
    print "blast.xml file exists, checked."
else:
    print "blast.xml file doesn't exist. Script will be terminated."
    sys.exit()

#BOTH FILES EXIST
#extract best blast hit info
command = "egrep -A 6 \"Iteration_query-def\" blast.xml | egrep \"Iteration_query-def|Hit_def\" | sed 's/Iteration_query-def//' | sed 's/Iteration_query-def//' | sed 's/Hit_def//'g | sed 's/<>//' | sed 's/<\/>//g' >blast.xml_parsed.txt"
process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True) 
process.communicate()[0] #Launch shell command

command = "bioawk -cfastx '{print $name \" \" gsub(s,s) \" \" length($seq) \" \" ((gsub(s,s))/(length($seq))*100)}' s='X' " + RMfile + " | sort -n | cut -d' ' -f1,4 >RM_parsed.txt"
process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True) 
process.communicate()[0] #Launch shell command


