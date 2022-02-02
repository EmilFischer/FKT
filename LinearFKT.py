import math
import numpy as np
import PlanarEmbedding
import networkx as nx
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

#--------- Create T1 ---------
t1 = planar.getSpanningTree()
t1Edges = list(t1.edges)

