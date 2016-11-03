#!python2
# Author: Kevin Sayers
# Assignment 5 implement a fasta parser and calculate the GC
# content of the sequences, determining the sequence with the
# highest GC content. The ID for the highest GC content is then
# output along with the GC content

# parse_fasta_file:
# input: a fasta file
# return: a dictionary containing the ID as the key and the sequence
# as the associated value. 

import sys #This import is just for command line argument parsing

def parse_fasta_file(fastafile):
	filelist = []
	sequences = {}

	with open(fastafile) as f:
		filelist = [line.strip() for line in f]

	i = 0
	# each time a > char is reached it becomes the current id 
	# it is assumed the following lines are valid nucleotide
	# sequences
	currentid = ''
	while i <= len(filelist)-1:
			if ">" in filelist[i]:
				sequences[filelist[i]] = ''
				currentid = filelist[i]
				i+=1
			else:
				sequences[currentid] = sequences[currentid] + filelist[i]
				i+=1

	return sequences

# get_gc_content: 
# input: a string of a  DNA nucleotide sequence
# return: the percentage of G/C in the sequence
def get_gc_content(sequence):
	gccount = 0.0 #a float so the percentage works later
	for nt in sequence:
		if nt == "G" or nt == "C":
			gccount += 1
	return gccount/len(sequence)*100




def main():
	# these if/else determines if the file was input as a 
	# command line argument or whether it should prompt
	if len(sys.argv) > 1:
		fastadata = parse_fasta_file(sys.argv[1])
	else:
		filein = raw_input("FASTA filename: ")
		fastadata = parse_fasta_file(filein)
	maxid = None
	maxgc = 0.0
	for i in fastadata:
		if get_gc_content(fastadata[i]) > maxgc:
			maxgc = get_gc_content(fastadata[i])
			maxid = i

	print maxid
	print maxgc

if __name__ == '__main__':
	main()
