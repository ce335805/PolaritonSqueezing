import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors
import h5py
from matplotlib import gridspec

import plotPolaritonFreqs

fontsize = 10

mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['font.size'] = 8  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.major.size'] = 3
mpl.rcParams['ytick.major.size'] = 3
mpl.rcParams['xtick.major.width'] = .7
mpl.rcParams['ytick.major.width'] = .7
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 8
mpl.rc('text', usetex=True)

mpl.rcParams['text.latex.preamble'] = [
    #    r'\renewcommand{\familydefault}{\sfdefault}',
    #    r'\usepackage[scaled=1]{helvet}',
    r'\usepackage[helvet]{sfmath}',
    #    r'\everymath={\sf}'
]


def plotTimeEvol2Bands():
    print("plotting more beautiful time evolution")

    # read in stuff
    file = h5py.File("../data/gsProp2BandsUa100Ub100Uud100Uss100epsA0epsB10gE20WD60FD1000NB6TS20WPh5000.hdf5", 'r')
    file2 = h5py.File("../data/gsProp2BandsUa100Ub100Uud100Uss100epsA0epsB10gE20WD70FD1000NB6TS20WPh5000.hdf5", 'r')

    #readInPrmsAndAssert(fileC, fileQ1)
    times = file['times'][()]
    pump = file['pump'][()]
    dOcc0 = file['dOcc0'][()]
    dOcc1 = file['dOcc1'][()]
    dOccUpDn = file['dOccUpDn'][()]
    dOccSigSig = file['dOccSigSig'][()]
    nPh = file['Nph1'][()]

    times2 = file2['times'][()]
    pump2 = file2['pump'][()]
    dOcc02 = file2['dOcc0'][()]
    dOcc12 = file2['dOcc1'][()]
    dOccUpDn2 = file2['dOccUpDn'][()]
    dOccSigSig2 = file2['dOccSigSig'][()]
    nPh2 = file2['Nph1'][()]


    dOccIntra = dOcc0 + dOcc1
    dOccInter = dOccUpDn + dOccSigSig
    dOccTot = dOcc0 + dOcc1 + dOccUpDn + dOccSigSig


    dOccIntra2 = dOcc02 + dOcc12
    dOccInter2 = dOccUpDn2 + dOccSigSig2
    dOccTot2 = dOcc02 + dOcc12 + dOccUpDn2 + dOccSigSig2


    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3., 2.)

    #ax.plot(times, pump)
    #ax.plot(times, nPh, linewidth = 0.5)
    #ax.plot(times, nPh2, linewidth = 0.5)
    ax.plot(times2, dOccTot, linewidth = .5)
    ax.plot(times2, dOccTot2, linewidth = .5)

    plt.savefig('./savedPlots/timeEvol.png', format='png', bbox_inches='tight', dpi=600)




