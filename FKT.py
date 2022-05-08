from decimal import Decimal, getcontext
import math
import PlanarEmbedding
import networkx as nx
import matplotlib.pyplot as plt
import NestedDissection
import Sparsification
import scipy.sparse as sparse
import MultiplyByTranspose as MbT
import os
import scipy.linalg
import time
from sympy import *
from scipy.linalg import lu
from decimal import *

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np

#--------- Create G ---------
G = nx.grid_2d_graph(20, 20)
A = nx.to_numpy_array(G)

planar = PlanarEmbedding.Planar(A)
G = planar.planar
Asparse = nx.adjacency_matrix(G)
A = nx.to_numpy_array(G)

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

facesList = []

while len(q) > 0:
  f = q.pop(0)

  if (f in used):
    continue
  
  facesList.append(f)

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

  if u != -1 and v != -1: #u and v are only -1 when f = outer face
    if (clockWiseEdges % 2 == 0):
      A[u, v] = 1
      A[v, u] = -1
    else:
      A[u, v] = -1
      A[v, u] = 1
  
  used.add(f)
  t2Leafs = set()
  del t2[f]
  for k in t2.keys():
    if f in t2[k]:
      t2[k].remove(f)
      if len(t2[k]) == 1:
        t2Leafs.add(k)
  q.extend(t2Leafs)

#--------- Output ---------
precision = G.number_of_nodes()
A = N(Matrix(A), precision)
#Naive computation of determinant
print("SYMPY DET COMPUTATION:")
start = time.time()
det = N(A.det(), precision)
end = time.time()

print("Elapsed time:", round(end - start, 3), "seconds")
print("Determinant:", int(det))
print("# of perf matches:", int(N(sqrt(det), precision)))
print("____________________________\n")

#Nested dissection computation of determinant
print("NESTED DISSECTION:")
start = time.time()
B = Sparsification.sparsify(A, Asparse)
BBT = B * B.transpose()
adjMatrix = np.array(BBT).astype(np.float64)
Gprime = nx.from_numpy_matrix(adjMatrix)

nd = NestedDissection.NestedDissection()
det = nd.determinant(Gprime, BBT, precision)
end = time.time()

print("Elapsed time:", round(end - start, 3), "seconds")
print("Determinant:", int(det))
print("# of perf matches:", int(sqrt(det)))

print("____________________________\n")

#n = np.shape(A)[0]
#for i in range(n):
#  for j in range(n):
#    if (A[i, j]) == 1:
#      A[j, i] = 0
#G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
#nx.draw_planar(G, with_labels=True)
#plt.show()