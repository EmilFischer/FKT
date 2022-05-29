from tkinter.ttk import Separator
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import PlanarEmbedding
import math
import sympy
from decimal import *
from scipy.linalg import lu
from sympy import *
import bisect

class NestedDissection:
    def __init__(self) -> None:
        self.numbers = dict()
        self.numbersInv = dict()
        self.fills = dict()

    def determinant(self, G, M, prec):
        self.number(G)
        self.fill(G)
        U = self.decomposition(M, prec)

        det = Decimal(1)
        for i in range(M.shape[0]):
            det *= U[(i, i)]
        det = abs(det)

        return det.sqrt()

    def separate(self, G) -> list:
        nodes = list(G.nodes)
        n = len(nodes)
        
        largestCC = max(nx.connected_components(G), key=len)
        source = next(iter(largestCC))

        levels = {source: 0}
        prevLvl = 0
        q = list(G.neighbors(source))
        parent = {source: source}
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

            for nbr in G.neighbors(vtx):
                if nbr not in parent:
                    parent[nbr] = vtx
                    q.append(nbr)

        A = []
        B = []
        S = []
        keys = levels.keys()
        for x in nodes:
            if x in self.numbers: continue
            if x in keys:
                if levels[x] == prevLvl:
                    S.append(x)
                else:
                    A.append(x)
            else:
                B.append(x)

        SNeighborhood = set()
        for v in S:
            for u in G.neighbors(v):
                if u in A:
                    SNeighborhood.add(u)
                    A.remove(u)
                elif u in B:
                    SNeighborhood.add(u)
                    B.remove(u)

        S.extend(SNeighborhood)

        return [A, B, S]

    def number(self, G):
        values = set()
        n = G.number_of_nodes()
        stack = [(G, 0, n-1)]
        used = set()

        alpha = 2./3
        beta = 2.83
        n0 = round(math.pow(beta/(1-alpha), 2))

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
                    if v not in self.numbers:
                        self.numbers[v] = a
                        self.numbersInv[a] = v
                        values.add(a)
                        a += 1
                        
            else:
                sets = self.separate(Gprime)
                i = len(sets[0])
                j = len(sets[1])
                k = len(sets[2])

                n = b-k+1
                Gprime = nx.Graph(Gprime)
                for v in sets[2]:
                    if v not in self.numbers: 
                        self.numbers[v] = n
                        self.numbersInv[n] = v
                        values.add(n)
                        n += 1
                        
                        nbrs = list(Gprime.neighbors(v))
                        for u in nbrs:
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

    def fill(self, G):
        n = G.number_of_nodes()
        for i in range(n):
            v = self.numbers[i]
            self.fills[v] = []
            for nbr in G.neighbors(i):
                if nbr > i:
                    nbr = self.numbers[nbr]
                    self.fills[v].append(nbr)

        for i in range(n):
            v = self.numbers[i]
            nbrs = self.fills[v]
            if len(nbrs) < 1:
                continue

            m = self.numbersInv[nbrs[0]]
            for vtx in nbrs:
                m = min(self.numbersInv[vtx], m)
            m = self.numbers[m]

            for w in nbrs:
                if w != m and w not in self.fills[m]:
                    self.fills[m].append(w)

    def decomposition(self, M, prec):
        fills = dict()
        for i in self.fills:
            if len(self.fills[i]) == 0: continue
            I = self.numbersInv[i]
            fills[I] = []
            for j in self.fills[i]:
                bisect.insort(fills[I], self.numbersInv[j])
        
        P = dict()
        for i in range(M.shape[0]):
            P[(i, i)] = Decimal(str(M[i, i]))
            if i in fills:
                for j in fills[i]:
                    P[(i, j)] = Decimal(str(M[i, j]))

        getcontext().prec = prec
        for i in fills:
            f = len(fills[i])
            for J in range(f):
                j = fills[i][J]
                s = P[(i, j)] / P[(i, i)]

                for L in range(J, f):
                    l = fills[i][L]
                    P[(j, l)] = P[(j, l)] - s * P[(i, l)]
            
            P[(i, j)] = s

        return P
