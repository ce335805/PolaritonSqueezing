import numpy as np
ID = np.array([[1., 0.], [0., 1.]])
JW = np.array([[-1., 0.], [0, 1.]])

cDag = np.array([[0., 1.], [0., 0.]])
cOp = np.array([[0., 0.], [1., 0.]])

N = 8
U = 1.
tHop = 1.

transposelist = np.zeros(2 * N, dtype = 'int')
