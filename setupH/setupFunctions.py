import numpy as np

from globals import *

def setupOpAtSite(site, localOp):

    C = np.zeros(0)

    for ind in np.arange(N):
        if (site == 0 and ind == 0):
            C = localOp
            continue
        if(ind < site):
            if(ind == 0):
                C = JW
            else:
                C = np.tensordot(JW, C, 0)
                #C = np.tensordot(C, JW, 0)
        elif(ind > site):
            C = np.tensordot(ID, C, 0)
            #C = np.tensordot(C, ID, 0)
        else:
            C = np.tensordot(localOp, C, 0)
            #C = np.tensordot(C, localOp, 0)

    C = np.transpose(C, transposelist)
    C = np.reshape(C, (2**N, 2**N))
    return C

def setupTransposeList():
    for ind in range(N):
        transposelist[ind] = int(2 * ind)
    for ind in range(N):
        transposelist[ind + N] = int(2 * ind + 1)