import numpy as np

from globals import *
from setupFunctions import *

eps0 = 0.
eps1 = 1.
tHop = -1.
U0 = 3.
U1 = 3.

def setupHTwoOrb():
    assert N == 8, "wrong number of sites for this Hamiltonian"

    HKin = setupHKin()
    HPot = onSitePotential()
    HInt = setupHInt()

    return HKin + HPot + HInt

def setupHKin():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2**N, 2**N), dtype = 'complex')

    for site in np.arange(N // 2):
        ci = setupOpAtSite(site, cOp)
        cip4 = setupOpAtSite(site + 4, cOp)
        ciDag = setupOpAtSite(site, cDag)
        cip4Dag = setupOpAtSite(site + 4, cDag)
        H += np.matmul(cip4Dag, ci) + np.matmul(ciDag, cip4)

    return tHop * H

def onSitePotential():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    for site in np.arange(N):
        ci = setupOpAtSite(site, cOp)
        ciDag = setupOpAtSite(site, cDag)
        if(site in [0, 1, 4, 5]):
            H += eps0 * np.matmul(ciDag, ci)
        else:
            H += eps1 * np.matmul(ciDag, ci)

    return

def onSitePotential0():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    for site in np.arange(N):
        ci = setupOpAtSite(site, cOp)
        ciDag = setupOpAtSite(site, cDag)
        if(site in [0, 1, 4, 5]):
            H += np.matmul(ciDag, ci)

    return H

def onSitePotential1():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    for site in np.arange(N):
        ci = setupOpAtSite(site, cOp)
        ciDag = setupOpAtSite(site, cDag)
        if(site in [2, 3, 6, 7]):
            H += np.matmul(ciDag, ci)

    return H


def setupHInt():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    #IDBig = np.identity(2 ** N)

    for site in np.arange(N // 2):
        ciUp = setupOpAtSite(2 * site, cOp)
        ciUpDag = setupOpAtSite(2 * site, cDag)
        ciDn = setupOpAtSite(2 * site + 1, cOp)
        ciDnDag = setupOpAtSite(2 * site + 1, cDag)

        if(site in [0, 2]):
            H += U0 *  np.matmul(np.matmul(ciUpDag, ciUp), np.matmul(ciDnDag, ciDn))
        else:
            H += U1 * np.matmul(np.matmul(ciUpDag, ciUp), np.matmul(ciDnDag, ciDn))

    return H


def interOrbitalHop0():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    for site in np.arange(2):
        ci = setupOpAtSite(site, cOp)
        cip2 = setupOpAtSite(site + 2, cOp)
        ciDag = setupOpAtSite(site, cDag)
        cip2Dag = setupOpAtSite(site + 2, cDag)
        H += np.matmul(cip2Dag, ci) + np.matmul(ciDag, cip2)

    return H

def interOrbitalHop1():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    for site in (np.arange(2) + 4):
        ci = setupOpAtSite(site, cOp)
        cip2 = setupOpAtSite(site + 2, cOp)
        ciDag = setupOpAtSite(site, cDag)
        cip2Dag = setupOpAtSite(site + 2, cDag)
        H += np.matmul(cip2Dag, ci) + np.matmul(ciDag, cip2)

    return H

#Docc in lower band
def setupDocc0():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    site = 0
    ciUp = setupOpAtSite(site, cOp)
    ciUpDag = setupOpAtSite(site, cDag)
    ciDn = setupOpAtSite(site + 1, cOp)
    ciDnDag = setupOpAtSite(site + 1, cDag)
    H += np.matmul(np.matmul(ciUpDag, ciUp), np.matmul(ciDnDag, ciDn))

    site = 4
    ciUp = setupOpAtSite(site, cOp)
    ciUpDag = setupOpAtSite(site, cDag)
    ciDn = setupOpAtSite(site + 1, cOp)
    ciDnDag = setupOpAtSite(site + 1, cDag)
    H += np.matmul(np.matmul(ciUpDag, ciUp), np.matmul(ciDnDag, ciDn))

    return H
#DOcc in upper band
def setupDocc1():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    site = 2
    ciUp = setupOpAtSite(site, cOp)
    ciUpDag = setupOpAtSite(site, cDag)
    ciDn = setupOpAtSite(site + 1, cOp)
    ciDnDag = setupOpAtSite(site + 1, cDag)
    H += np.matmul(np.matmul(ciUpDag, ciUp), np.matmul(ciDnDag, ciDn))

    site = 6
    ciUp = setupOpAtSite(site, cOp)
    ciUpDag = setupOpAtSite(site, cDag)
    ciDn = setupOpAtSite(site + 1, cOp)
    ciDnDag = setupOpAtSite(site + 1, cDag)
    H += np.matmul(np.matmul(ciUpDag, ciUp), np.matmul(ciDnDag, ciDn))

    return H

def interOrbitalIntUpDn():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    for site in [0, 4]:
        ciUp = setupOpAtSite(site, cOp)
        ciUpDag = setupOpAtSite(site, cDag)
        ciDn = setupOpAtSite(site + 3, cOp)
        ciDnDag = setupOpAtSite(site + 3, cDag)
        H += np.matmul(np.matmul(ciUpDag, ciUp), np.matmul(ciDnDag, ciDn))

    for site in [1, 5]:
        ciDn = setupOpAtSite(site, cOp)
        ciDnDag = setupOpAtSite(site, cDag)
        ciUp = setupOpAtSite(site + 1, cOp)
        ciUpDag = setupOpAtSite(site + 1, cDag)
        H += np.matmul(np.matmul(ciUpDag, ciUp), np.matmul(ciDnDag, ciDn))

    return H


def interOrbitalIntSigSig():
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    for site in [0, 1, 4, 5]:
        ciUp0 = setupOpAtSite(site, cOp)
        ciUp0Dag = setupOpAtSite(site, cDag)
        ciUp1 = setupOpAtSite(site + 2, cOp)
        ciUp1Dag = setupOpAtSite(site + 2, cDag)
        H += np.matmul(np.matmul(ciUp0Dag, ciUp0), np.matmul(ciUp1Dag, ciUp1))

    return H

def ni(site):
    assert N == 8, "wrong number of sites for this Hamiltonian"
    H = np.zeros((2 ** N, 2 ** N), dtype='complex')

    cSite = setupOpAtSite(site, cOp)
    cSiteDag = setupOpAtSite(site, cDag)
    H += np.matmul(cSiteDag, cSite)

    return H

def nc0():
    return ni(0) + ni(1)

def nd0():
    return ni(2) + ni(3)

def nc1():
    return ni(4) + ni(5)

def nd1():
    return ni(6) + ni(7)



