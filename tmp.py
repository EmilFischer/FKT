import numpy as np
a = np.array([[1,2], [3,4]])
zeros = np.zeros((8, 8))
z=zeros.shape[0]
n=a.shape[0]

zeros[:-z+n,:-z+n] = a
print(zeros)