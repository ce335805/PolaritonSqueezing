import numpy as np
import h5py

ID = np.array([[1, 0], [0, 1]])
JW = np.array([[1, 0], [0, -1]])

cOp = np.array([[0, 0], [1, 0]])
cDag = np.array([[0, 1], [0, 0]])
N = 4
U = 1.
tHop = -1.

transposelist = np.zeros(2 * N, dtype = 'int')

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
        elif(ind > site):
            C = np.tensordot(ID, C, 0)
        else:
            C = np.tensordot(localOp, C, 0)

    C = np.transpose(C, transposelist)
    C = np.reshape(C, (2**N, 2**N))
    return C

def setupTransposeList():
    for ind in range(N):
        transposelist[ind] = int(2 * ind)
    for ind in range(N):
        transposelist[ind + N] = int(2 * ind + 1)


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
        H += np.matmul(cN, c0) + np.matmul(c0Dag, cN)

    return tHop * H

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

    return tHop * H

def main():
    print("Let's setup H!")

    setupTransposeList()

    print("transposeList = {}".format(transposelist))

    #op = setupOpAtSite(0, ID)
    #print(op)

    HHop = setUpHoppingH()
    HInd = setUpHInt()
    H = HInd + HHop
    #print(H)

    filename = "HN{}U{}.hdf5".format(N, int(U * 10))
    print("H.shape = {}".format(H.shape))

    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data = np.real(H), dtype = 'double')
    f.create_dataset("Imag", data = np.imag(H), dtype = 'double')
    f.close()


    dOcc = setUpDOcc()
    print("dOcc.shape = {}".format(dOcc.shape))

    filename = "dOccN{}.hdf5".format(N)

    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data = np.real(dOcc), dtype = 'double')
    f.create_dataset("Imag", data = np.imag(dOcc), dtype = 'double')
    f.close()


    couplingH = setUpCouplingH()
    print("couplingH.shape = {}".format(couplingH.shape))

    filename = "couplingN{}.hdf5".format(N)

    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data = np.real(couplingH), dtype = 'double')
    f.create_dataset("Imag", data = np.imag(couplingH), dtype = 'double')
    f.close()

main()