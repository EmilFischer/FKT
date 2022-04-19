from decimal import Decimal
import math
import numpy as np
import PlanarEmbedding
import networkx as nx
import matplotlib.pyplot as plt
import NestedDissection
import Sparsification
import scipy.sparse as sparse
import MultiplyByTranspose as MbT
from scipy.linalg import lu

#--------- Create G ---------
G = nx.grid_2d_graph(12, 12)
Asparse = nx.adjacency_matrix(G)
A = Asparse.todense()
planar = PlanarEmbedding.Planar(A)

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
#planar = PlanarEmbedding.Planar(A)
#G = planar.G
#Asparse = nx.adjacency_matrix(G)

#--------- Create T1 ---------
t1 = planar.getSpanningTree()
t1Edges = list(t1.edges)

#--------- Find edges that are in G but not in T1 ---------
#--------- while also orienting T1                ---------
notInT1 = set()
for v1, v2 in planar.edges:
  if (v1, v2) not in t1Edges and (v2, v1) not in t1Edges:
    notInT1.add((v1, v2))
    notInT1.add((v2, v1))
    # Set edges that are not in T1 to 0 in A (orient later)
    A[v1, v2] = 0
    A[v2, v1] = 0
  else:
    # Orient all edges in T1
    A[v2, v1] = -1
    A[v1, v2] = 1

#--------- Create T2 ---------
t2 = dict()
q = [planar.detectFaces()]
used = set()
while len(q) > 0:
  f = q.pop(0)
  if (f in used):
    continue
  
  t2[f] = set()
  edges = f.getEdges()
  for adj in f.getAdjFaces():
    for e in adj.getEdges():
      if e in notInT1 and e in edges: 
        t2[f].add(adj)
  
  used.add(f)
  q.extend(f.getAdjFaces())

#--------- Find leafs of T2 ---------
t2Leafs = set()
for k in t2.keys():
  if len(t2[k]) == 1:
    t2Leafs.add(k)

#--------- Loop through T2, pfaffian orienting the graph ---------
q = list(t2Leafs)
used = set()
while len(q) > 0:
  f = q.pop(0)
  if (f in used):
    continue

  vts = f.getVertices()
  clockWiseEdges = 0
  u = -1
  v = -1
  for i in range(len(vts)):
    v1 = vts[i-1]
    v2 = vts[i]
    if A[v1, v2] == 1:
      clockWiseEdges += 1
    elif A[v1, v2] == 0:
      u = v1
      v = v2

  if u != -1 and v != -1:
    if (clockWiseEdges % 2 == 0):
      A[u, v] = 1
      A[v, u] = -1
    else:
      A[u, v] = -1
      A[v, u] = 1
  
  used.add(f)
  q.extend(f.getAdjFaces())

#--------- Output ---------
#Naive computation of determinant
det = round(np.linalg.det(A))
print("NAIVE APPROACH:")
print("Determinant:", det)
print("# of perf matches:", round(math.sqrt(det)))
print("____________________________\n")

#Nested dissection computation of determinant
print("NESTED DISSECTION:")
B = Sparsification.sparsify(A, Asparse)
BBt = np.matmul(B, B.transpose())

Gprime = nx.from_numpy_matrix(BBt)
nd = NestedDissection.NestedDissection()
det = nd.determinant(Gprime, BBt)

print("Determinant:", det)
print("# of perf matches:", round(math.sqrt(det)))
print()

#n = np.shape(A)[0]
#for i in range(n):
#  for j in range(n):
#    if (A[i, j]) == 1:
#      A[j, i] = 0
#G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
#nx.draw_planar(G, with_labels=True)
#plt.show()