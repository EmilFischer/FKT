import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import PlanarEmbedding
import math

class NestedDissection:
    T1 = None

    def __init__(self) -> None:
        self.numbers = dict()
        self.fills = dict()

    def determinant(self, G):
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

    treeIter = 0
    lastparent = None
    parent = None
    fillLeft = True
    def createTreeNode(self, A, B, separator):
        node = Node(A, B, separator)

        if self.treeIter == 0:
            self.T1 = node
            self.parent = self.T1
            self.lastParent = self.T1
        elif self.treeIter % 2 == 1:
            self.parent.child1 = node
        else:
            self.parent.child2 = node
            tmp = self.parent
            if self.fillLeft:
                self.parent = self.lastParent.child1
            else:
                self.parent = self.lastParent.child2
                
            self.lastParent = tmp
            self.fillLeft = not self.fillLeft

        self.treeIter += 1

    def number(self, G):
        alpha = 2./3
        beta = 6
        n0 = math.pow(beta/(1-alpha), 2)
        values = set()

        q = [(G, 0, G.number_of_nodes()-1)]
        used = set()
        while len(q) > 0:
            triple = q.pop(0)
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

                self.createTreeNode(sets[0], sets[1], sets[2])

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
                q.append((sub, a, a+i-1))

                union = sets[2]
                union.extend(sets[1])
                sub = Gprime.subgraph(union)
                q.append((sub, b-k-j+1, b-k))
    
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
        M = nx.adjacency_matrix(G)
        LU = np.zeros((n, n))
        for i in self.fills.keys():
            for j in range(len(self.fills[i])):
                k = self.fills[i][j]
                s = M[i, k]
                
                for l in range(j, len(self.fills[i])):
                    LU[k, l] = M[j, l] - s * M[i, l]
                
                LU[i, k] = s
        return LU
            
class Node:
    child1 = None
    child2 = None
    A = None
    B = None
    separator = None

    def __init__(self, A, B, separator):
        self.A = A
        self.B = B
        self.separator = separator

def sparseRepresentation(M):
    n = np.shape(M)[0]
    r = []
    for i in range(n):
        r.append([])
        for j in range(n):
            if M[i, j] != 0:
                r[i].append((j, M[i, j]))
    return r

#with open('100x100grid.npy', 'rb') as f:
#    A = np.load(f)

G = nx.grid_2d_graph(2, 2)
M = nx.to_scipy_sparse_matrix(G)

nd = NestedDissection()
dissection = nd.determinant(G)