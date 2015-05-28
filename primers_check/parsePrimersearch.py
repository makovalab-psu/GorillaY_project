#parsePrimersearch.py primersearch_file file_with_primers
import sys
import re

primersearch_file=sys.argv[1]
primers_file=sys.argv[2]
primers_order=[]
primers=dict()
forward=dict()
reverse=dict()
amplicons=[]

with open(primers_file, "r") as file:
	for file_line in file:
				line=file_line.rstrip()
				arr=line.split('\t')
				#print arr
				contig_name=arr[0]
				F=arr[1]
				R=arr[2]
				forward[contig_name]=F
				reverse[contig_name]=R
				primers_order.append(contig_name)


with open(primersearch_file, "r") as file:
			next_key=""
			for file_line in file:
				line=file_line.rstrip()
				#print line
				if ("Primer name" in line):
					arr=line.split(" ")
					key=next_key
					next_key=arr[2]
					#print ("key: " + key)

					if (key!=""):
						amplicons.sort()
						primers[key]=amplicons
						amplicons=[]
				else:
					if ("Amplimer length" in line):
						arr=line.split(" ")
						value=int(arr[2])
						#print ("value: " + value)
						amplicons.append(value)
			if line is not None:
				arr=line.split(" ")
				key=next_key
				next_key=arr[2]
				#print ("key: " + key)

				if (key!=""):
					amplicons.sort()
					primers[key]=amplicons
					amplicons=[]


#print primers
print ("Contig_name\tsmallest\tlessThan20kbp\tforward\treverse\tnumberOfAmplicons\tfirst_10_values")
for key in primers_order:
	values=primers[key]
	smallest_formatted="NaN"
	smallest=-1
	if (len(values)>0):
		smallest=int(min(values))
		smallest_formatted=re.sub("(\d)(?=(\d{3})+(?!\d))", r"\1,", "%d" % smallest)
	star="X"
	if (smallest<20000 and smallest>0):
		star="*"
	print (key + "\t" + str(smallest_formatted) + "\t" + star + "\t" + forward[key] + "\t" + reverse[key] + "\t" + str(len(values)) + "\t" + str(values[:10]))



