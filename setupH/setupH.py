import numpy as np
import h5py

from globals import *
from setupFunctions import *
import singleBandHubbard
import twoOrbitalH

def projectList(val, op):

    deleteInd = set()
    for diagInd, diagEntry in enumerate(np.diag(op)):
        #print(diagEntry)
        if(diagEntry != val):
            deleteInd.add(diagInd)

    return deleteInd

def setupN():
    nOp = np.zeros((2**N, 2**N), dtype = 'complex')
    for site in np.arange(N):
        ci = setupOpAtSite(site, cOp)
        ciDag = setupOpAtSite(site, cDag)
        nOp += np.matmul(ciDag, ci)
    return nOp

def setupNupMNdn():
    op = np.zeros((2**N, 2**N), dtype = 'complex')
    for site in np.arange(N//2):
        ciUp = setupOpAtSite(2 * site, cOp)
        ciDagUp = setupOpAtSite(2 * site, cDag)
        ciDn = setupOpAtSite(2 * site + 1, cOp)
        ciDagDn = setupOpAtSite(2 * site + 1, cDag)
        op += np.matmul(ciDagUp, ciUp) - np.matmul(ciDagDn, ciDn)
    return op


def projectN(N, H):
    nOp = setupN()
    nUpMNDn = setupNupMNdn()

    deleteIndN = projectList(N, nOp)
    deleteIndNUpMNDn = projectList(0, nUpMNDn)
    deleteIndTot = deleteIndN.union(deleteIndNUpMNDn)
    deleteIndArr = np.array(list(deleteIndTot), dtype='int')

    Hcut = np.delete(H, deleteIndArr, axis = 0)
    Hcut = np.delete(Hcut, deleteIndArr, axis = 1)

    return Hcut

def writeH2Orb(Hcut, interOrb0, interOrb1, dOcc0, dOcc1):
    filename = "./savedOperators/HU0_{}_U1_{}_eps0_{}_eps1_{}.hdf5".format(twoOrbitalH.U0, twoOrbitalH.U1, twoOrbitalH.eps0, twoOrbitalH.eps1)

    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(Hcut), dtype='double')
    f.create_dataset("Imag", data=np.imag(Hcut), dtype='double')
    f.close()

    filename = "./savedOperators/InterOrb0.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(interOrb0), dtype='double')
    f.create_dataset("Imag", data=np.imag(interOrb0), dtype='double')
    f.close()

    filename = "./savedOperators/InterOrb1.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(interOrb1), dtype='double')
    f.create_dataset("Imag", data=np.imag(interOrb1), dtype='double')
    f.close()


    filename = "./savedOperators/dOcc0.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(dOcc0), dtype='double')
    f.create_dataset("Imag", data=np.imag(dOcc0), dtype='double')
    f.close()


    filename = "./savedOperators/dOcc1.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(dOcc1), dtype='double')
    f.create_dataset("Imag", data=np.imag(dOcc1), dtype='double')
    f.close()

def main():
    print("Let's setup H!")

    setupTransposeList()

    print("transposeList = {}".format(transposelist))

    #setupOpAtSite(1, cDag)
    #exit()

    #mat1 = np.array([[1, 2], [3, 4]])
    #mat2 = np.array([[5, 6], [7, 8]])
    #ten1 = np.tensordot(mat1, mat2, 0)
    #np.transpose(ten1, transposelist)
    #ten1 = np.reshape(ten1, (4, 4))
    #print(ten1)
    #print()
    ##print(np.tensordot(mat2, mat1, 0))
    #exit()

    #HHop = singleBandHubbard.setUpHoppingH()
    #HInd = singleBandHubbard.setUpHInt()
    #H = HInd + HHop
    ##print(np.real(H))
    #Hcut = projectN(2, H)
    #print(np.real(Hcut))
    #exit()

    particleN = 2

    H = twoOrbitalH.setupHTwoOrb()
    Hcut = projectN(particleN, H)
    interOrb0 = twoOrbitalH.interOrbitalHop0()
    interOrb0Cut = projectN(particleN, interOrb0)
    interOrb1 = twoOrbitalH.interOrbitalHop1()
    interOrb1Cut = projectN(particleN, interOrb1)

    dOcc0 = twoOrbitalH.setupDocc0()
    dOcc0Cut = projectN(particleN, dOcc0)
    dOcc1 = twoOrbitalH.setupDocc1()
    dOcc1Cut = projectN(particleN, dOcc1)

    print(np.real(np.diag(Hcut)))

    eVals, _ = np.linalg.eigh(Hcut)
    print(eVals)

    writeH2Orb(Hcut, interOrb0Cut, interOrb1Cut, dOcc0Cut, dOcc1Cut)
    exit()

    filename = "HN{}U{}.hdf5".format(N, int(U * 10))
    print("H.shape = {}".format(H.shape))

    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data = np.real(Hcut), dtype = 'double')
    f.create_dataset("Imag", data = np.imag(Hcut), dtype = 'double')
    f.close()


    #dOcc = setUpDOcc()
    #print("dOcc.shape = {}".format(dOcc.shape))
#
    #filename = "dOccN{}.hdf5".format(N)
#
    #f = h5py.File(filename, 'w')
    #f.create_dataset("Real", data = np.real(dOcc), dtype = 'double')
    #f.create_dataset("Imag", data = np.imag(dOcc), dtype = 'double')
    #f.close()
#
#
    #couplingH = setUpCouplingH()
    #print("couplingH.shape = {}".format(couplingH.shape))
#
    #filename = "couplingN{}.hdf5".format(N)
#
    #f = h5py.File(filename, 'w')
    #f.create_dataset("Real", data = np.real(couplingH), dtype = 'double')
    #f.create_dataset("Imag", data = np.imag(couplingH), dtype = 'double')
    #f.close()

main()