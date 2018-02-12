import numpy as np
from timeit import default_timer as timer



def numpy_read(N, A):

    for i in range(0,N):
        b = A[55,55]

def list_read(N, B):

    for i in range(0,N):
        b = A[55][55]

N = 1000000
A = np.zeros((1000,1000))
B = [[0] * 1000] * 1000

start = timer()
numpy_read(N, A)
end = timer()
print("numpy_read time: %f" % (end-start))

start = timer()
list_read(N, B)
end = timer()
print("list_read time: %f" % (end-start))