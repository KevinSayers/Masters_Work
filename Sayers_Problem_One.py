#Kevin Sayers
# 1/7/2017
from Bio import SeqIO
import itertools

#read_fasta: reads the fasta file, a tuple is created containing each sequence and fasta id. 
#inputs: fasta file path
#outputs: none
#return: list of tuples with sequence and id
def read_fasta(infile):
	return [(str(i.seq), str(i.id)) for i in SeqIO.parse(infile,"fasta")]

#check_overlap: checks if the suffic of length k matches the prefix of length k returns true if so.
#inputs: two strings, and the k length integer
#outputs: none
#return: boolean
def check_overlap(v,w,k):
	return v[-k:] == w[0:k]

#find_edges: iterates over all combinations of fasta ids and finds those that overlap.
#inputs: list of tuples with sequence and fasta id, the length integer k
#outputs: none
#return: the edge list
def find_edges(fastas, k):
	edge_list = []
	for i in fastas:
		for j in fastas:
			if i[1] != j[1]:
				if check_overlap(i[0],j[0],k):
					edge_list.append(i[1] + " " + j[1])

	return edge_list


def main():
	k = 3
	recs = read_fasta("r1.txt")
	results = find_edges(recs,k)

	outfile = open("results.txt","w")
	for i in results:
		print (i)
		outfile.write(i + "\n")

	outfile.close()


if __name__ == "__main__":
	main()