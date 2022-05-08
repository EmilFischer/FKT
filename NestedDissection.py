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
        n = shape(M)[0]
        print("n:", n)
        print("Edges:", G.number_of_edges())

        self.number(G)
        print("Numbers found!")
        print("Length of numbers:", len(self.numbers))
        c = 0
        for i in range(n):
            if self.numbers[i] != i:
                c += 1
        print("Different numbers:", c)

        self.fill(G)
        print("Fill-ins found!")

        c = 0
        q = 0
        for i in range(n):
            l = len(self.fills[i])
            c += l
            q = max(q, l)
        print("Edges + fill-ins:", c)
        print("Most edges + fill-ins in row:", q)

        U = self.decomposition(M, prec)
        print("Decomposition found!")

        #_,U2,_ = M.LUdecomposition()
        #U2 = N(U2, prec)

        det = Float(1.0, prec)
        for i in range(U.shape[0]):
            #if U2[i, i] != U[i, i]: print(i, N(U2[i, i], 10), N(U[i, i], 10))
            det = N(det * U[i, i], prec)
        det = round(det)

        return N(sqrt(Abs(det)), prec)

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
        stack = [(G, 0, G.number_of_nodes()-1)]
        used = set()

        alpha = 2./3
        beta = 6
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
                        while a in values:
                            a += 1
                        self.numbers[v] = a
                        self.numbersInv[a] = v
                        values.add(a)
            else:
                sets = self.separate(Gprime)
                i = len(sets[0])
                j = len(sets[1])
                k = len(sets[2])

                n = b-k+1
                Gprime = nx.Graph(Gprime)
                for v in sets[2]:
                    if v not in self.numbers:
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

    def decomposition(self, M, precision):
        fills = dict()
        for i in self.fills:
            if len(self.fills[i]) == 0: continue
            I = self.numbersInv[i]
            fills[I] = []
            for j in self.fills[i]:
                bisect.insort(fills[I], self.numbersInv[j])
        
        for i in range(len(fills)):
            for j in fills[i]:
                s = N(M[i, j] / M[i, i], precision)

                for l in fills[i]:
                    if l < j: continue
                    M[j, l] = N(M[j, l] - s * M[i, l], precision)
            
            M[i, j] = s

        return N(M, precision)
