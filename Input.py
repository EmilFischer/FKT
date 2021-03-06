import networkx as nx
import numpy as np
from FKT import FKT
from Halin import halin

print("___________________________________")
print("_____________INITIATED_____________")
print("___________________________________\n")
# Set graph size
for n in range(8, 16, 2):
    print("n:", n)
    #--------- Create G ---------
    #-Halin graph-
    #G = halin(3, n)

    #-Grid graph-
    G = nx.grid_2d_graph(n, n)

    #-Path Graph-
    #G = nx.path_graph(n)

    #--------- Run FKT ---------
    #Optional parameters to state which 
    #methods to use to compute determinant:
    #nesteddissection = True/False
    #naive = True/False
    #sympy = True/False
    
    FKT(G, naive=True)


#-Graph from adj. matrix-
#A = np.array([
#  [0,1,1,0,1,0,0,1],
#  [1,0,0,1,1,1,0,0],
#  [1,0,0,1,0,0,1,1],
#  [0,1,1,0,0,1,1,0],
#  [1,1,0,0,0,0,0,0],
#  [0,1,0,1,0,0,0,0],
#  [0,0,1,1,0,0,0,0],
#  [1,0,1,0,0,0,0,0]
#])
#G = nx.from_numpy_matrix(A)
#FKT(G)