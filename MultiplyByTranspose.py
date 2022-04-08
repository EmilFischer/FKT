import scipy.sparse as sparse
# https://stackoverflow.com/questions/24566633/which-is-the-best-way-to-multiply-a-large-and-sparse-matrix-with-its-transpose

def MultiplyByTranspose(A):
    I = []
    J = []
    V = []
    n = A.shape[0]
    for i in range(n):
        for j in range(i, n):
            X = 0.0
            for k in range(n):
                X += A[i][k] * A[k][j]
            if X != 0.0:
                I.append(i)
                J.append(j)
                V.append(X)
                I.append(j)
                J.append(i)
                V.append(X)
    
    return sparse.coo_matrix((V,(I,J)),shape=(n, n))