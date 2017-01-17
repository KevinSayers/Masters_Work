#Kevin Sayers
#1/15/2017
from Bio.Seq import Seq
import networkx as nx #used to plot the network graph 
import matplotlib.pyplot as plt

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

#getcyclic: finds the cyclic superstring. It finds the edges needed for each kmer once only.
#this is then passed to recurseedges which traverses the nodes creating the superstring.
#input: edges, kmers, and k value
#output: none
#return: a shortest cyclic superstring 
def getcyclic(edges,kmers,k):
	visited = []
	sups = []
	superstr = ""

	finaledges = []
	for i in kmers:
		if (i[:k-1],i[1:]) in edges:
			finaledges.append((i[:k-1],i[1:]))
	return recurseedges(finaledges)

#recurseedges: this function traverses the linked nodes creating the shortest superstring.
#input: edges
#output: none
#return: a shortest cyclic superstring
def recurseedges(edges):
	curi = 0
	superstr = ""
	#The algorithm works by removing visited edges as it traverses. 
	#it is initialized at an index of 0 
	while len(edges) > 0:
		current = edges.pop(curi)
		superstr += current[1][-1]
		for i in range(0,len(edges)):
			if edges[i][0] == current[1]:
				curi = i
	return (superstr)

#plotgraph: This function uses networkx module to plot the nodes and edges. 
#input: nodes, edges, k
#output: a matplotlib graph of the network
#return: none
def plotgraph(nodes,edges,k):
	finaledges = []
	for i in nodes:
		if (i[:k-1],i[1:]) in edges:
			finaledges.append((i[:k-1],i[1:]))
	newnodes = []
	finaledges = sorted(finaledges)
	edgelabels = []
	labeldict = {}

	for p in finaledges:
		newnodes.append(p[0])
		labeldict[p] = p[1][-1]

	g = nx.DiGraph()
	g.add_nodes_from(newnodes)
	g.add_edges_from(finaledges)
	pos = nx.circular_layout(g)
	# labels = nx.draw_networkx_edge_labels(g,pos,labeldict)
	nx.draw(g, with_labels=True)
	plt.show()

				


def main():
	input = readfile("n4.txt")
	k = len(input[0])
	elements = getelements(input)
	result = getedges(elements,k)
	final = getcyclic(result,input,k)
	out = open("q4.txt","w")
	out.write(final)
	plotgraph(input,result,k)
	# result = sorted(result)
	# for i in result:
	# 	print (i)
	# 	out.write(str(i).replace("'","")+"\n")
	# out.close()


if __name__ == "__main__":
	main()