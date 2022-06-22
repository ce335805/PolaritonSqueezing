import numpy as np

from globals import *
from setupFunctions import *

def setUpHoppingH():

    H = np.zeros((2**N, 2**N), dtype = 'complex')

    for site in np.arange(N - 2):
        ci = setupOpAtSite(site, cOp)
        cip2 = setupOpAtSite(site + 2, cOp)
        ciDag = setupOpAtSite(site, cDag)
        cip2Dag = setupOpAtSite(site + 2, cDag)
        H += np.matmul(cip2Dag, ci) + np.matmul(ciDag, cip2)

    if(N > 4):
        cN = setupOpAtSite(N - 2, cOp)
        c0 = setupOpAtSite(0, cOp)
        cNDag = setupOpAtSite(N - 2, cDag)
        c0Dag = setupOpAtSite(0, cDag)
        H += np.matmul(cNDag, c0) + np.matmul(c0Dag, cN)

        cN = setupOpAtSite(N - 1, cOp)
        c0 = setupOpAtSite(1, cOp)
        cNDag = setupOpAtSite(N - 1, cDag)
        c0Dag = setupOpAtSite(1, cDag)
        H += np.matmul(cNDag, c0) + np.matmul(c0Dag, cN)

    return H

def setUpHInt():
    H = np.zeros((2**N, 2**N), dtype = 'complex')

    for site in np.arange(N // 2):
        ciUp = setupOpAtSite(2 * site, cOp)
        ciUpDag = setupOpAtSite(2 * site, cDag)
        ciDown = setupOpAtSite(2 * site + 1, cOp)
        ciDownDag = setupOpAtSite(2 * site + 1, cDag)

        IDBig = np.identity(2**N)

        H += np.matmul(np.matmul(ciUpDag, ciUp) - 0.5 * IDBig, np.matmul(ciDownDag, ciDown) - 0.5 * IDBig)

    return U * H

def setUpDOcc():
    dOcc = np.zeros((2**N, 2**N), dtype = 'complex')

    for site in np.arange(N // 2):
        ciUp = setupOpAtSite(2 * site, cOp)
        ciUpDag = setupOpAtSite(2 * site, cDag)
        ciDown = setupOpAtSite(2 * site + 1, cOp)
        ciDownDag = setupOpAtSite(2 * site + 1, cDag)

        IDBig = np.identity(2**N)

        dOcc += np.matmul(np.matmul(ciUpDag, ciUp) - 0.5 * IDBig, np.matmul(ciDownDag, ciDown) - 0.5 * IDBig)

    return dOcc


def setUpCouplingH():

    H = np.zeros((2**N, 2**N), dtype = 'complex')

    for site in np.arange(N - 2):
        ci = setupOpAtSite(site, cOp)
        cip2 = setupOpAtSite(site + 2, cOp)
        ciDag = setupOpAtSite(site, cDag)
        cip2Dag = setupOpAtSite(site + 2, cDag)
        H += np.matmul(cip2Dag, ci) - np.matmul(ciDag, cip2)

    if(N > 4):
        cN = setupOpAtSite(N - 2, cOp)
        c0 = setupOpAtSite(0, cOp)
        cNDag = setupOpAtSite(N - 2, cDag)
        c0Dag = setupOpAtSite(0, cDag)
        H += np.matmul(cNDag, c0) - np.matmul(c0Dag, cN)

        cN = setupOpAtSite(N - 1, cOp)
        c0 = setupOpAtSite(1, cOp)
        cNDag = setupOpAtSite(N - 1, cDag)
        c0Dag = setupOpAtSite(1, cDag)
        H += np.matmul(cNDag, c0) - np.matmul(c0Dag, cN)

    return H