import math
import numpy as np
import galois
import networkx as nx
import nxmetis
import PlanarEmbedding
import matplotlib.pyplot as plt
from sympy import *

#--------- Create G ---------
A = np.array([
  [0,1,1,0,1,0,0,0,0,0,0,1],
  [1,0,0,1,0,1,1,0,0,0,0,0],
  [1,0,0,1,0,0,0,0,0,1,1,0],
  [0,1,1,0,0,0,0,1,1,0,0,0],
  [1,0,0,0,0,1,0,0,0,0,0,1],
  [0,1,0,0,1,0,1,0,0,0,0,0],
  [0,1,0,0,0,1,0,1,0,0,0,0],
  [0,0,0,1,0,0,1,0,1,0,0,0],
  [0,0,0,1,0,0,0,1,0,1,0,0],
  [0,0,1,0,0,0,0,0,1,0,1,0],
  [0,0,1,0,0,0,0,0,0,1,0,1],
  [1,0,0,0,1,0,0,0,0,0,1,0]
])
planar = PlanarEmbedding.Planar(A)

#--------- Find all faces in the graph ---------
planar.detectFaces()
faces = planar.facesList

#--------- Create orientation matrix ---------
edgeIdx = dict()
idx = -1
def getEdgeIdx(e):
  if e in edgeIdx.keys():
    return edgeIdx[(e[0], e[1])]
  else:
    global idx
    idx += 1
    edgeIdx[(e[0], e[1])] = idx
    return idx

width = int(len(planar.edges)/2)
height = len(faces)
GF = galois.GF(2)
orientationMatrix = GF.Zeros((height, width))

rhs = GF.Zeros(height)

i = 0
for f in faces:
  u = -1
  excess = 1
  for v in range(len(f)):
    vtx1 = f.vertices[u]
    vtx2 = f.vertices[v]
    if vtx1 < vtx2:
      j = getEdgeIdx((vtx1, vtx2))
      orientationMatrix[i, j] = 1
    else: 
      excess -= 1
      j = getEdgeIdx((vtx2, vtx1))
      orientationMatrix[i, j] = 1
    u += 1
  rhs[i] = excess % 2
  i += 1 

print("Orientation Matrix:\n", orientationMatrix)

#result = np.linalg.solve(orientationMatrix, rhs)
#print(result)

#--------- Output ---------
#G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
G = nx.from_numpy_matrix(A)
nx.draw_planar(G, with_labels=True)
plt.show()
