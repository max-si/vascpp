
import numpy as np
from scipy.sparse import csr_matrix

# nnz = 20
# numVessels = 7
# indptr = np.array([0, 2, 4, 8, 11, 14, 17, 20])
# print(len(indptr))
# indices = np.array([2, 8, 7, 13, 0, 3, 4, 8, 2, 4, 5, 2, 3, 6, 3, 6, 7, 4, 5, 7])
# print(len(indices))
# data = np.ones(nnz)

# A = csr_matrix((data, indices, indptr), shape=(numVessels, numVessels)).toarray()

# print(A)


# ParMETIS Example

# Serial CSR
nnz = 44
numNodes = 15
indptr = np.array([0, 2, 5, 8, 11, 13, 16, 20, 24, 28, 31, 33, 36, 39, 42, 44])
indices = np.array([1, 5, 0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9, 0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14, 5, 11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13])
print(len(indptr))
print(len(indices))
data = np.ones(nnz)
A = csr_matrix((data, indices, indptr), shape=(numNodes, numNodes)).toarray()
print(A)

# Distributed CSR Processor 0
nnz = 13
numVessels = 5
indptr = np.array([0,2, 5, 8, 11, 13])
indices = np.array([1,5,0,2,6,1,3,7,2,4,8,3,9])
data = np.ones(nnz)
B = csr_matrix((data, indices, indptr), shape=(numVessels, numVessels)).toarray()
print(B)

# Distributed CSR Processor 1
nnz = 13
numVessels = 5
indptr = np.array([0,2, 5, 8, 11, 13])
indices = np.array([1,5,0,2,6,1,3,7,2,4,8,3,9])
data = np.ones(nnz)
C = csr_matrix((data, indices, indptr), shape=(numVessels, numVessels)).toarray()
print(C)