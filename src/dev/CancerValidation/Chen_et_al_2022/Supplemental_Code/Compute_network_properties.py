# Functions to calculate network properties

# Construct a protein-protein interaction network
import networkx as nx
edges = set()
# Define a set of high-quality human protein-protein interactions (the HINT database)
hint = set([tuple(sorted(line.strip("\n").split("\t")[:2])) for line in open("HomoSapiens_binary_hq.txt").readlines()[1:]])
# Intersect the HINT interactome with interface predictions to define the 3D interacome network
for line in open("H_sapiens_interfacesALL.txt"):
    p1,p2=line.strip().split("\t")[:2]
    if tuple(sorted([p1,p2])) not in hint: continue
    edges.add(tuple(sorted([p1,p2])))   
G = nx.Graph(list(edges))

# Pre-compute betweenness values for network G
node_betweenness = nx.betweenness_centrality(G)
edge_betweenness = nx.edge_betweenness_centrality(G)

# Compute network degree of a given protein in G
def networkDegree(prot):
    return len(G[prot])

# Compute node betweenness of a given protein in G
def nodeBetweenness(prot):
    return node_betweenness[prot]

# Compute edge betweenness of a given pair of proteins (edge) in G
def edgeBetweenness([p1,p2]):
    return edge_betweenness[tuple(sorted([p1,p2]))]

