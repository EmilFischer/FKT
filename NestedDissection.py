import numbers
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import PlanarEmbedding
import math

class NestedDissection:
    def __init__(self, G) -> None:
        self.numbers = dict()
        self.number(G)
        
        self.fills = dict()
        self.fill(G)

        self.decomposition(G)

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

            if lvl > prevLvl and len(levels) >= n/2:
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

        return [A, B, S]

    def number(self, G):
        alpha = 2./3
        beta = 6
        n0 = math.pow(beta/(1-alpha), 2)
        values = set()

        stack = [(G, 0, G.number_of_nodes()-1)]
        used = set()
        while len(stack) > 0:
            triple = stack.pop()
            Gprime = triple[0]
            nodes = list(Gprime.nodes)

            graphHash = str(nodes)
            if graphHash in used: continue
            used.add(graphHash)

            a = triple[1]
            b = triple[2]

            if len(nodes) <= n0:
                for v in nodes:
                    if v not in self.numbers.keys():
                        while a in values:
                            a += 1
                        self.numbers[v] = a
                        values.add(a)
                        
            else:
                sets = self.separate(Gprime)
                i = len(sets[0])
                j = len(sets[1])
                k = len(sets[2])

                iter = 0
                for n in range(b-k+1, b):
                    v = sets[2][iter]
                    if v not in self.numbers.keys():
                        if n not in values:
                            self.numbers[v] = n
                            values.add(n)
                            continue
                        for u in Gprime.neighbors(v):
                            if u in sets[2]:
                                Gprime.remove_edge(u, v)
                    iter += 1

                union = sets[2].copy()
                union.extend(sets[0])
                sub = Gprime.subgraph(union)
                stack.append((sub, a, a+i-1))

                union = sets[2]
                union.extend(sets[1])
                sub = Gprime.subgraph(union)
                stack.append((sub, b-k-j+1, b-k))
    
    def fill(self, G):
        nodes = list(G.nodes)
        for v in self.numbers.values():
            nbrs = list(G.neighbors(nodes[v]))
            m = self.numbers[nbrs[0]]
            for vtx in nbrs:
                m = min(self.numbers[vtx], m)

            if m not in self.fills:
                self.fills[m] = []
            
            for vtx in nbrs:
                w = self.numbers[vtx]
                if w != m:
                    self.fills[m].append(w)
            
    def decomposition(self, G):
        n = G.number_of_nodes()
        L = np.zeros((n, n))
        print(self.fills)
        for i in self.fills.keys():
            for j in self.fills[i]:
                s = j/(float(i)+1)
                for l in range(len(self.fills[i])):
                    L[j, l] = L[j, l] - s * L[i, l]
                L[i, j] = s
        print(L)
            
#with open('100x100grid.npy', 'rb') as f:
#    A = np.load(f)

#planar = PlanarEmbedding.Planar(A)
G = nx.grid_2d_graph(2, 2)
dissection = NestedDissection(G)