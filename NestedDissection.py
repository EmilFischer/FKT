from tkinter.ttk import Separator
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import PlanarEmbedding
import math
from decimal import *
from scipy.linalg import lu
import bisect

class NestedDissection:
    T1 = None

    def __init__(self) -> None:
        self.numbers = dict()
        self.numbersInv = dict()
        self.fills = dict()

    def determinant(self, G, M):
        _, _, U = lu(M)
        actualFillIns = []
        for i in range(M.shape[0]):
            for j in range(i, M.shape[0]):
                if U[i, j] != 0 and M[i, j] == 0:
                    actualFillIns.append((i, j))


        self.number(G)

        for i in self.numbers:
            if i < 0 or i > 343:
                print("Oh no!", i)

        self.fill(G)
        
        c = 0
        for p in actualFillIns:
            if p[1] not in self.fills[p[0]]:
                c += 1
        print("# unidentified fill ins:", c)

        U = self.decomposition(M)

        det = Decimal(1)
        for i in range(U.shape[0]):
            det *= Decimal(U[i, i])
        
        return math.isqrt(round(abs(det)))

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

        SNbrs = set()
        for v in S:
            for nbr in G.neighbors(v):
                if nbr not in S:
                    SNbrs.add(nbr)
                    if nbr in A:
                        A.remove(nbr)
                    elif nbr in B:
                        B.remove(nbr)
        S.extend(SNbrs)

        return [A, B, S]

    def number(self, G):
        alpha = 2./3
        beta = 6
        n0 = round(math.pow(beta/(1-alpha), 2))
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
                        self.numbersInv[a] = v
                        values.add(a)
                        
            else:
                Gprime = nx.Graph(Gprime)
                sets = self.separate(Gprime)
                i = len(sets[0])
                j = len(sets[1])
                k = len(sets[2])

                n = b-k+1
                for v in sets[2]:
                    if v not in self.numbers.keys():
                        while n in values:
                            n += 1
                            
                        self.numbers[v] = n
                        self.numbersInv[n] = v
                        values.add(n)

                        nbrs = list(Gprime.neighbors(v))
                        for u in nbrs:
                            if u in sets[2]:
                                Gprime.remove_edge(u, v)

                union = sets[2].copy()
                union.extend(sets[1])
                sub = Gprime.subgraph(union)
                stack.append((sub, b-k-j+1, b-k))

                union = sets[2].copy()
                union.extend(sets[0])
                sub = Gprime.subgraph(union)
                stack.append((sub, a, a+i-1))
    
    def fill(self, G):
        n = G.number_of_nodes()
        for i in range(n):
            self.fills[i] = []
            for nbr in G.neighbors(i):
                if nbr > i:
                    self.fills[i].append(nbr)
        
        for i in range(n):
            v = self.numbers[i]
            nbrs = self.fills[v]
            if len(nbrs) < 1:
                continue

            m = nbrs[0]
            for vtx in nbrs:
                m = min(vtx, m)

            for w in nbrs:
                if w != m and w not in self.fills[m]:
                    bisect.insort(self.fills[m], w)
            
    def decomposition(self, M):
        for i in range(M.shape[0]):
            for j in self.fills[i]:
                s = M[i, j] / M[i, i]

                for l in set(self.fills[i] + self.fills[j]):
                    if l < j: continue
                    M[j, l] = M[j, l] - s * M[i, l]
                
                M[i, j] = s
        
        return M
