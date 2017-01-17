#Kevin Sayers
# 1/7/2017
from Bio import SeqIO
import itertools

#read_fasta: returns a list of the sequences from a fasta file
#inputs: fasta filename
#outputs: none
#return: list of sequences from file
def read_fasta(infile):
	return [str(i.seq) for i in SeqIO.parse(infile,"fasta")]

#check_overlap: returns the overlap between two strings
#inputs: two strings
#outputs: none
#return: returns either the overlap or if no overlap 0
def check_overlap(v,w): 
	for i in range(len(v),0,-1):
		if v[-i:] == w[:i]:
			return i
	return 0

#construct_superstring: recursively combines pairs of strings that match
#by at least 50%. Once only one string is remaining, this superstring is 
#returned.
#inputs: list of sequences
#outputs: none
#return: either a recursive call with the a new list of fastas after merging
#or the superstring 
def construct_superstring(fastas):
	results = []
	for j in range(0,len(fastas)):
		max_overlap = 0
		max_pos = None
		temp = fastas[j]
		for i in range(0,len(fastas)):
			if fastas[i] != temp:
				overlap = check_overlap(temp,fastas[i])
				if overlap > max_overlap:
					max_overlap = overlap
					max_pos = i
		if max_overlap > len(fastas[j]) * 0.5:
			results.append(temp + str(fastas[max_pos][max_overlap:]))

	if len(results) > 1:
		return construct_superstring(results)
	else:
		return results[0]


def main():
	recs = read_fasta("r2.txt")
	results = construct_superstring(recs)
	print (results)
	outfile = open("out2.txt","w")
	outfile.write(results)

if __name__ == "__main__":
	main()