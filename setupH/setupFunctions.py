import numpy as np

from globals import *

def setupOpAtSite(site, localOp):

    opArr = np.zeros((N, 2, 2), dtype = complex)

    for ind in np.arange(N):
        if(ind < site):
            opArr[ind, :, :] = JW
        elif(ind > site):
            opArr[ind, :, :] = ID
        else:
            opArr[ind, :, :] = localOp

    retOp = np.tensordot(opArr[0, :, :], opArr[1, :, :], 0)
    for ind in np.arange(N-2):
        retOp = np.tensordot(retOp, opArr[ind + 2, :, :], 0)

    #print(np.real(opArr))
    #print()
    #print(np.real(retOp))

    retOp = np.transpose(retOp, transposelist)
    retOp = np.reshape(retOp, (2**N, 2**N))
    #print(np.real(retOp))
    return retOp


#def setupOpAtSite(site, localOp):
#
#    C = np.zeros(0)
#
#    for ind in np.arange(N):
#        if (site == 0 and ind == 0):
#            C = localOp
#            continue
#        if(ind < site):
#            if(ind == 0):
#                C = JW
#            else:
#                #C = np.tensordot(JW, C, 0)
#                C = np.tensordot(C, JW, 0)
#        elif(ind > site):
#            #C = np.tensordot(ID, C, 0)
#            C = np.tensordot(C, ID, 0)
#        else:
#            #C = np.tensordot(localOp, C, 0)
#            C = np.tensordot(C, localOp, 0)
#
#    C = np.transpose(C, transposelist)
#    C = np.reshape(C, (2**N, 2**N))
#    return C

def setupTransposeList():
    for ind in range(N):
        transposelist[ind] = int(2 * ind)
    for ind in range(N):
        transposelist[ind + N] = int(2 * ind + 1)