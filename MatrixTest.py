import numpy as np
import math
import scipy.linalg
from scipy.linalg import lu
from decimal import Decimal, getcontext
from sympy import Matrix

#Random 100x100 matrix with integers from 1 to 4
with open('test.npy', 'rb') as f:
  A = np.load(f)

print("NUMPY DET COMPUTATION:")
det = np.linalg.det(A)
det = abs(int(det))

print("Determinant:", det)
print("____________________________\n")

print("SCIPY DET COMPUTATION:")
det = scipy.linalg.det(A)
det = abs(int(det))

print("Determinant:", det)
print("____________________________\n")

print("SCIPY LU DECOMPOSITION:")
_, _, U = lu(A)
getcontext().prec = 100
det = Decimal(1)
for i in range(U.shape[0]):
  det *= Decimal(U[i, i])
det = abs(det)

print("Determinant:", int(det))
print("____________________________\n")

print("SYMPY DET COMPUTATION")
A = Matrix(A)
det = abs(A.det())

print("Determinant:", int(det))
print("____________________________\n")