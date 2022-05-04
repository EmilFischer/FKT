from decimal import Decimal, getcontext
import math
import PlanarEmbedding
import networkx as nx
import matplotlib.pyplot as plt
import NestedDissection
import Sparsification
import scipy.sparse as sparse
import MultiplyByTranspose as MbT
import os
import scipy.linalg
from sympy import *
from scipy.linalg import lu
from decimal import *

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np

n = 64
T = np.zeros((n, n), dtype=np.longdouble)

for j in range(n):
    for k in range(n):
        a = math.pow(4 * (math.cos(pi * (j+1) / (n+1))), 2)
        b = math.pow(4 * (math.cos(pi * (k+1) / (n+1))), 2)
        T[j, k] = math.pow((a + b), 0.25)

R = np.prod(np.prod(T))
print(R)