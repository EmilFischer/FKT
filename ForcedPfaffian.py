import math
import numpy as np
import networkx as nx
import PlanarEmbedding
import matplotlib.pyplot as plt

#--------- Create G ---------
A = np.array([
  [0, 1, 0, 1, 0, 0],
  [1, 0, 1, 1, 0, 0],
  [0, 1, 0, 0, 1, 1],
  [1, 1, 0, 0, 1, 0],
  [0, 0, 1, 1, 0, 1],
  [0, 0, 1, 0, 1, 0.]
])
planar = PlanarEmbedding.Planar(A)

#--------- Set lower triangular to all zeros ---------
n = np.shape(A)[0]
for i in range(n):
  for j in range(i):
    A[i, j] = 0

#--------- Loop through faces and orient ---------
q = [planar.detectFaces()]
for f in planar.facesList: print(f)
used = set()
while len(q) > 0:
  f = q.pop(0)
  used.add(f)

  sum = 0
  e = (0, 0)
  i = -1
  vts = f.vertices
  for j in range(len(vts)):
    u = vts[i]
    v = vts[j]
    if A[u, v] == 1:
      sum += 1
      e = (u, v)
    i += 1

  checkAdj = False
  if sum % 2 == 0:
    checkAdj = True
    A[e[0], e[1]] = 0
    A[e[1], e[0]] = 1

  if checkAdj:
    q.extend(f.adjFaces)
  else:
    for adj in f.adjFaces:
      if adj not in used:
        q.append(adj)

#--------- Output ---------
G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
nx.draw_planar(G, with_labels=True)
plt.show()
