import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

class face:
  def __init__(self):
    self.verticies = []
    self.adjFaces = set()

  def __str__(self) -> str:
    return self.verticies.__str__()

  def addVerticies(self, vts):
    self.verticies = vts

  def getVerticies(self) -> list:
    return self.verticies
  
  def addAdjFace(self, f):
    self.adjFaces.add(f)

  def getAdjFaces(self) -> set:
    if (self in self.adjFaces):
      self.adjFaces.remove(self)
    return self.adjFaces

  def getEdges(self) -> set:
    edges = set()
    n = len(self.verticies)
    if n < 2: return edges
    for i in range(1, n):
      j = i-1
      edges.add((self.verticies[i], self.verticies[j]))
      edges.add((self.verticies[j], self.verticies[i]))
    
    edges.add((self.verticies[-1], self.verticies[0]))
    edges.add((self.verticies[0], self.verticies[-1]))
    return edges


class Planar:
  def __init__(self, matrix):
    self.G = nx.from_numpy_matrix(matrix)
    self.n = np.shape(matrix)[0]
    self.planar = nx.check_planarity(self.G)[1]
    self.edges = list(self.planar.edges)

  def draw(self):
    nx.draw_planar(self.planar, with_labels=True)
    plt.show()

  def getSpanningTree(self):
    return nx.minimum_spanning_tree(self.planar)

  def detectFaces(self):
    nbrs = self.planar.get_data()
    f = face()
    faces = dict()

    # Ignore all single edges/vertices
    singleVts = set()
    for v in nbrs:
      if (len(nbrs[v]) == 1):
        singleVts.add(v)
        next = nbrs[v][0]
        
        nextNbrs = nbrs[next]
        nextNbrs = [ vtx for vtx in nextNbrs if vtx not in singleVts ]

        while (len(nextNbrs) == 1):
          singleVts.add(next)
          next = nextNbrs[0]
          nextNbrs = nbrs[next]
          nextNbrs = [ vtx for vtx in nextNbrs if vtx not in singleVts ]
          
    # Error handling for graph with no faces
    if len(singleVts) == self.n:
      f.addVerticies = list(singleVts)
      return f

    # Detect each face by looping over each vertex
    for origin in range(self.n):
      if (origin in singleVts): 
        continue

      idx = 0
      currNbrs = nbrs[origin]
      next = currNbrs[idx]
      while ((origin, next) in faces.keys()) or (next in singleVts):
        idx -= 1
        next = currNbrs[idx]
        if (idx == -len(currNbrs)):
          break
      if ((origin, next) in faces.keys()):
        continue

      curr = origin
      verticies = [origin]
      f = face()

      while (next != origin):
        verticies.append(next)

        if ((next, curr) in faces.keys()):
          adj = faces[(next, curr)]
          adj.addAdjFace(f)
          f.addAdjFace(adj)
        
        faces[(curr, next)] = f
        
        idx = nbrs[next].index(curr)-1
        curr = next
        next = nbrs[next][idx]
        while next in singleVts:
          idx -= 1
          next = nbrs[curr][idx]

      if (next, curr) in faces.keys():
        adj = faces[(next, curr)]
        adj.addAdjFace(f)
        f.addAdjFace(adj)

      faces[(curr, next)] = f
      n = len(verticies)
      if (n > 2): 
        f.addVerticies(verticies)

    return f
