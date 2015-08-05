'''
arguments are R1_postrecY_reads, R2_reference, R2_postrecY_reads_to_build

'''
import sys


word_count_reference = 100

pattern_file = open(sys.argv[1],'r')
reference_file = open(sys.argv[2],'r')
complement_file = open(sys.argv[3],'w')

found = 0
curr_ref_line = []*1
big_pat_line = []*1
saved_record = []*4
small_pat_line = []*1
small_curr_ref_line = []*1
do = 1
tries = 0


for big_pat_line in pattern_file : 
	if big_pat_line[0:4] == '@HIS' :  
		small_pat_line = big_pat_line.split(' ')[0]
		if found  == 1 :
			#write_saved_record
			complement_file.write(saved_record[0])
			complement_file.write(saved_record[1])
			complement_file.write(saved_record[2])
			complement_file.write(saved_record[3])
			saved_record = []*4
			found = 0
		while found == 0 :
			big_curr_ref_line = reference_file.readline()
			small_curr_ref_line = big_curr_ref_line.split(' ')[0]
			#print curr_ref_line
			#print curr_ref_line[0]
			if small_curr_ref_line[0:4] == '@HIS' :
			#if do == 1 :
				#print small_pat_line
				#print small_curr_ref_line
				
				if small_pat_line == small_curr_ref_line :
					#print small_pat_line
				 	#print small_curr_ref_line
					#save 4 linesi
					#print "match, yo!"
					saved_record.append(big_curr_ref_line)
					saved_record.append(reference_file.readline())
					saved_record.append(reference_file.readline())
					saved_record.append(reference_file.readline())
					#print saved_record
					is_first_time = 0
					found = 1
					tries+=1
				else : 
					ignore=reference_file.readline()
					ignore=reference_file.readline()
					ignore=reference_file.readline()
		



complement_file.write(saved_record[0])
complement_file.write(saved_record[1])
complement_file.write(saved_record[2])
complement_file.write(saved_record[3])


pattern_file.close()
reference_file.close()
complement_file.close()