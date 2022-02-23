import networkx as nx
import numpy as np
import PlanarEmbedding
import math

class NestedDissection:
    def __init__(self, G) -> None:
        self.numbers = dict()
        nodes = list(G.nodes)
        self.number(G, 1, len(nodes), 0)

    def separate(self, G) -> list:
        nodes = list(G.nodes)
        source = nodes[0]
        n = len(nodes)
        
        levels = {source: 0}
        prevLvl = 0
        q = list(G.neighbors(source))
        parent = dict()
        parent[source] = source
        for x in q:
            parent[x] = source

        while len(q) > 0:
            vtx = q.pop(0)

            p = parent[vtx]
            lvl = levels[p] + 1

            if lvl > prevLvl and len(levels) > n/2:
                break

            levels[vtx] = lvl
            prevLvl = lvl

            nbrs = list(G.neighbors(vtx))
            keys = parent.keys()
            for x in nbrs:
                if x not in keys:
                    parent[x] = vtx
                    q.append(x)

        A = []
        B = []
        S = []
        keys = levels.keys()
        for x in nodes:
            if x in keys:
                if levels[x] == prevLvl:
                    S.append(x)
                else:
                    A.append(x)
            else:
                B.append(x)

        print(prevLvl)
        print(levels)
        return [A, B, S]


    def number(self, G, a, b):
        alpha = 0
        beta = 0
        nodes = list(G.nodes)

        if len(nodes) <= math.pow(beta/(1-alpha), 2):
            for v in nodes:
                if v not in self.numbers.keys():
                    self.numbers[v] = a
                    a += 1
        else:
            sets = self.separate(G)
            i = len(sets[0])
            j = len(sets[1])
            k = len(sets[2])

            iter = 0
            for n in range(b-k+1, b):
                v = sets[2][iter]
                iter += 1
                if v not in self.numbers.keys():
                    self.numbers[v] = n
                    for u in G.Neighbors(v):
                        if u in sets[2]:
                            G.remove_edge(u, v)

            Gprime = G
            for v in nodes:
                if v in sets[0]:
                    Gprime.remove_node(v)
            self.number(Gprime, b-k-j+1, b-k)

            Gprime = G
            for v in nodes:
                if v in sets[1]:
                    Gprime.remove_node(v)
            self.number(Gprime, b-k-j-i+1, a+i-1)

            
    

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
dissection = NestedDissection(planar.G)
