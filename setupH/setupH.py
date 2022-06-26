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
    
def write2OrbOps(Hkin, onsitePot0, onsitePot1, interOrb0, interOrb1, dOcc0, dOcc1, intOrbUpDn, intOrbSigSig):
    filename = "./savedOperators/HKin.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(Hkin), dtype='double')
    f.create_dataset("Imag", data=np.imag(Hkin), dtype='double')
    f.close()

    filename = "./savedOperators/onsitePot0.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(onsitePot0), dtype='double')
    f.create_dataset("Imag", data=np.imag(onsitePot0), dtype='double')
    f.close()

    filename = "./savedOperators/onsitePot1.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(onsitePot1), dtype='double')
    f.create_dataset("Imag", data=np.imag(onsitePot1), dtype='double')
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

    filename = "./savedOperators/intOrbUpDn.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(intOrbUpDn), dtype='double')
    f.create_dataset("Imag", data=np.imag(intOrbUpDn), dtype='double')
    f.close()

    filename = "./savedOperators/intOrbSigSig.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(intOrbSigSig), dtype='double')
    f.create_dataset("Imag", data=np.imag(intOrbSigSig), dtype='double')
    f.close()

def main():
    setupTransposeList()
    print("transposeList = {}".format(transposelist))

    particleN = 2

    Hkin = twoOrbitalH.setupHKin()
    HkinCut = projectN(particleN, Hkin)
    onsitePot0 = twoOrbitalH.onSitePotential0()
    onsitePot0Cut = projectN(particleN, onsitePot0)
    onsitePot1 = twoOrbitalH.onSitePotential1()
    onsitePot1Cut = projectN(particleN, onsitePot1)
    interOrb0 = twoOrbitalH.interOrbitalHop0()
    interOrb0Cut = projectN(particleN, interOrb0)
    interOrb1 = twoOrbitalH.interOrbitalHop1()
    interOrb1Cut = projectN(particleN, interOrb1)

    dOcc0 = twoOrbitalH.setupDocc0()
    dOcc0Cut = projectN(particleN, dOcc0)
    dOcc1 = twoOrbitalH.setupDocc1()
    dOcc1Cut = projectN(particleN, dOcc1)

    intOrbUpDn = twoOrbitalH.interOrbitalIntUpDn()
    intOrbUpDnCut = projectN(particleN, intOrbUpDn)
    intOrbSigSig = twoOrbitalH.interOrbitalIntSigSig()
    intOrbSigSigCut = projectN(particleN, intOrbSigSig)

    write2OrbOps(HkinCut, onsitePot0Cut, onsitePot1Cut, interOrb0Cut, interOrb1Cut, dOcc0Cut, dOcc1Cut, intOrbUpDnCut, intOrbSigSigCut)

main()