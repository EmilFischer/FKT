import networkx as nx
import numpy as np
import PlanarEmbedding
import math

class NestedDissection:
    def __init__(self, G) -> None:
        self.numbers = dict()
        self.number(G)

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

        return [A, B, S]

    def number(self, G):
        alpha = 2./3
        beta = 8
        n0 = math.ceil(math.pow(beta/(1-alpha), 2))
        print(n0)

        stack = [(G, 1, G.number_of_nodes())]
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
                        self.numbers[v] = a
                        a += 1
            else:
                sets = self.separate(Gprime)
                i = len(sets[0])
                j = len(sets[1])
                k = len(sets[2])

                iter = 0
                for n in range(b-k+1, b):
                    v = sets[2][iter]
                    iter += 1
                    if v not in self.numbers.keys():
                        self.numbers[v] = n
                        for u in Gprime.neighbors(v):
                            if u in sets[2]:
                                Gprime.remove_edge(u, v)

                union = sets[2].copy()
                union.extend(sets[0])
                sub = Gprime.subgraph(union)
                stack.append((sub, a, a+i-1))

                union = sets[2]
                union.extend(sets[1])
                sub = Gprime.subgraph(union)
                stack.append((sub, b-k-j+1, b-k))

            
with open('100x100grid.npy', 'rb') as f:
    A = np.load(f)

planar = PlanarEmbedding.Planar(A)
dissection = NestedDissection(planar.G)

print("length:", len(dissection.numbers))
for n in dissection.numbers.values():
    if n < 0 or n > 9999:
        print(n)