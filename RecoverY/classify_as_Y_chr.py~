'''
arguments are R1_reference_reads kmer_table_after_threshold R1_postrecY_reads_to_build

'''

import sys

# User HAS TO set kmer_size and s1
kmer_size = 11 

# The higher you set s1, the more strict it is going to be on autosomes 
# but you will also lose some Y
# this is for 250bp reads

s1 = 125



input_all_types_reads = open(sys.argv[1],'r')
surely_male_kmers = open(sys.argv[2],'r')
classified_as_male_reads = open(sys.argv[3],'w')

surely_male_kmers_set = set (line[:kmer_size] for line in surely_male_kmers)

flag = 0
save_readname = ['']
save_read = ['']
save_rest = ['']
prospective_matches = 0

line_count = 0
record_count = 0
read_count = 0

for line in input_all_types_reads :
	line_count+=1
	if line[0:4] == '@HIS' :
		record_count +=1 
		if flag == 100 :
			classified_as_male_reads.write(str(save_readname))
			classified_as_male_reads.write(str(save_read))
			for element in save_rest : 
				classified_as_male_reads.write(str(element)+'\n')
			#print save_readname
			#print save_read
			#print save_rest
		flag = 1
		save_readname = line[:]
		save_read =['']
		save_rest = ['']

	elif line[0] != '@' and flag == 1 :
		read_count += 1
		start_index = 0
		prospective_matches = 0 
		while (start_index < (len(line)-25)) :
			current_kmer = line[start_index:start_index+kmer_size]
			#print current_kmer
			if current_kmer in surely_male_kmers_set :
				#print 'found a match!'
				#print line
				#print current_kmer
				prospective_matches += 1
			start_index += 1

		if prospective_matches > s1 :
			save_read = line[:-1]
			flag = 100
		else :
			flag = 0
	elif line[0:4] != '@' and flag==100 :
		save_rest.append(line[:-1])
	
	

#print line_count
#print record_count
#print read_count

input_all_types_reads.close()
surely_male_kmers.close()
classified_as_male_reads.close()
