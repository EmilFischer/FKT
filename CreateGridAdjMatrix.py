import numpy as np

def make_matrix(rows, cols):
    n = rows*cols
    M = np.zeros((n, n))
    for r in range(rows):
        for c in range(cols):
            i = r*cols + c
            # Two inner diagonals
            if c > 0: M[i-1, i] = M[i, i-1] = 1
            # Two outer diagonals
            if r > 0: M[i-cols, i] = M[i, i-cols] = 1
    return M

width = 100
height = 100

M = make_matrix(height, width)
fileName = str(height) + 'x' + str(width) + 'grid.npy'
with open(fileName, 'wb') as f:
    np.save(f, M)