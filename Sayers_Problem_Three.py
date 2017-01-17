#Kevin Sayers
#1/15/2017
from Bio.Seq import Seq

#readfile: reads data from file
#input: filename
#output: none
#return: list of lines from the file
def readfile(infile):
	data = []
	with open(infile) as f:
		data = f.read().splitlines()
	return data

#reverse_comp: returns the reverse complement of a sequence
#input: a DNA sequence as a string
#output: none
#return: a string representing the reverse complement
def reverse_comp(seq):
	return str(Seq(seq).reverse_complement())

#getelements: makes a set of kmers and their reverse complement
#input: list of kmers
#output: none
#return: a set of all kmers and their reverse complement
def getelements(inlist):
	result_set = set()
	for i in inlist:
		result_set.add(i)
		result_set.add(reverse_comp(i))
	return result_set

#getedges: returns the edges for each node
#input: kmers, and k value
#output:none
#return: a list of tuples for the edges
def getedges(elements,k):
	result = []
	for i in elements:
		result.append((i[0:k-1],i[1:k]))

	return result



def main():
	input = readfile("r3.txt")
	k = len(input[0])
	elements = getelements(input)
	result = getedges(elements,k)

	out = open("q3.txt","w")
	result = sorted(result)
	for i in result:
		print (i)
		out.write(str(i).replace("'","")+"\n")
	out.close()


if __name__ == "__main__":
	main()