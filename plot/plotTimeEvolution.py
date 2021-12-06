import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors
import h5py

fontsize = 14

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
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['figure.titlesize'] = 8
mpl.rc('text', usetex=True)

mpl.rcParams['text.latex.preamble'] = [
    #    r'\renewcommand{\familydefault}{\sfdefault}',
    #    r'\usepackage[scaled=1]{helvet}',
    r'\usepackage[helvet]{sfmath}',
    #    r'\everymath={\sf}'
]


def main():
    print("plotting some beautiful time evolution")

    #read in stuff
    fileC = h5py.File("../data/timeEvolutionResultsClassicalDrive.hdf5", 'r')
    fileQ = h5py.File("../data/timeEvolutionResultsQuantumDrive.hdf5", 'r')
    times = fileC['times'][()]
    pump = fileC['pump'][()]
    dOccC = fileC['dOcc'][()]
    XptC = fileC['Xpt'][()]
    XptSqrC = fileC['XptSqr'][()]
    NptC = fileC['Npt'][()]
    XphC = fileC['Xph'][()]
    XphSqrC = fileC['XphSqr'][()]
    NphC = fileC['Nph'][()]

    fileQ = h5py.File("../data/timeEvolutionResultsQuantumDrive.hdf5", 'r')
    dOccQ = fileQ['dOcc'][()]
    XptQ = fileQ['Xpt'][()]
    XptSqrQ = fileQ['XptSqr'][()]
    NptQ = fileQ['Npt'][()]
    XphQ = fileQ['Xph'][()]
    XphSqrQ = fileQ['XphSqr'][()]
    NphQ = fileQ['Nph'][()]

    print("times.shape = {}".format(times.shape))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    #ax.plot(times, pump, color = 'cornflowerblue' , label = 'pump')
    #ax.plot(times, NphC, color = 'olive' , label = r'$N_{\rm phon}$ - Classical')
    #ax.plot(times, NphQ, color = 'rosybrown' , label = r'$N_{\rm phon}$ - Quantum')
    #ax.plot(times, NptC, color = 'mediumseagreen' , label = r'$N_{\rm phot}$')
    #ax.plot(times, XphSqrC * 0.2, label = r'$\langle X_{\rm phon}^2 \rangle$', color = 'olive')
    #ax.plot(times, XptSqrC * 0.1, label = r'$\langle X_{\rm phot}^2 \rangle$', color = 'rosybrown')
    #ax.plot(times, XphC, label = r'$\langle X_{\rm phon} \rangle$', color = 'olive')
    #ax.plot(times, XptC, label = r'$\langle X_{\rm phot} \rangle$', color = 'rosybrown')

    ax.plot(times, dOccC + .5, color = 'olive', label = 'dOcc - Classical')
    #ax.plot(times, dOccQ + .5, color = 'rosybrown', label = 'dOcc - Quantum')

    plt.legend()

    #ax.set_ylim(-1e10, 1e10)

    plt.show()



main()
