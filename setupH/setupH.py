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
    
def write2OrbOps(Hkin, onsitePot0, onsitePot1, interOrb0, interOrb1, dOcc0, dOcc1, intOrbUpDn, intOrbSigSig, n0, n1, nc0, nd0, nc1, nd1):
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

    filename = "./savedOperators/n0.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(n0), dtype='double')
    f.create_dataset("Imag", data=np.imag(n0), dtype='double')
    f.close()

    filename = "./savedOperators/n1.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(n1), dtype='double')
    f.create_dataset("Imag", data=np.imag(n1), dtype='double')
    f.close()

    filename = "./savedOperators/nc0.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(nc0), dtype='double')
    f.create_dataset("Imag", data=np.imag(nc0), dtype='double')
    f.close()

    filename = "./savedOperators/nd0.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(nd0), dtype='double')
    f.create_dataset("Imag", data=np.imag(nd0), dtype='double')
    f.close()

    filename = "./savedOperators/nc1.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(nc1), dtype='double')
    f.create_dataset("Imag", data=np.imag(nc1), dtype='double')
    f.close()

    filename = "./savedOperators/nd1.hdf5"
    print("writing to file " + filename)
    f = h5py.File(filename, 'w')
    f.create_dataset("Real", data=np.real(nd1), dtype='double')
    f.create_dataset("Imag", data=np.imag(nd1), dtype='double')
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

    n0 = twoOrbitalH.ni(0)
    n0cut = projectN(particleN, n0)
    n1 = twoOrbitalH.ni(2)
    n1cut = projectN(particleN, n1)

    nc0 = twoOrbitalH.nc0()
    nc0cut = projectN(particleN, nc0)
    nd0 = twoOrbitalH.nd0()
    nd0cut = projectN(particleN, nd0)
    nc1 = twoOrbitalH.nc1()
    nc1cut = projectN(particleN, nc1)
    nd1 = twoOrbitalH.nd1()
    nd1cut = projectN(particleN, nd1)

    write2OrbOps(HkinCut,
                 onsitePot0Cut,
                 onsitePot1Cut,
                 interOrb0Cut,
                 interOrb1Cut,
                 dOcc0Cut,
                 dOcc1Cut,
                 intOrbUpDnCut,
                 intOrbSigSigCut,
                 n0cut,
                 n1cut,
                 nc0cut,
                 nd0cut,
                 nc1cut,
                 nd1cut)

main()