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


def plotGSProps():
    print("plotting GS props")

    fileOnePh = h5py.File("../data/gsProps1PhG20N6.hdf5", 'r')
    times = fileOnePh['times'][()]
    pump = fileOnePh['pump'][()]
    dOccOnePh = fileOnePh['dOcc'][()]
    XptOnePh = fileOnePh['Xpt'][()]
    XptSqrOnePh = fileOnePh['XptSqr'][()]
    NptOnePh = fileOnePh['Npt'][()]
    XphOnePh = fileOnePh['X1ph'][()]
    XphSqrOnePh = fileOnePh['X1phSqr'][()]
    NphOnePh = fileOnePh['N1ph'][()]

    fileTwoPh = h5py.File("../data/gsProps2PhG20N6.hdf5", 'r')
    times = fileTwoPh['times'][()]
    pump = fileTwoPh['pump'][()]
    dOccTwoPh = fileTwoPh['dOcc'][()]
    XptTwoPh = fileTwoPh['Xpt'][()]
    XptSqrTwoPh = fileTwoPh['XptSqr'][()]
    NptTwoPh = fileTwoPh['Npt'][()]
    XphTwoPh = fileTwoPh['X1ph'][()]
    XphSqrTwoPh = fileTwoPh['X1phSqr'][()]
    NphTwoPh = fileTwoPh['N1ph'][()]


    fig = plt.figure()
    ax = fig.add_subplot(111)


    ax.plot(times, dOccOnePh + 0.5, color='olive', label='dOcc, $\omega_{\rm P}$ - $X^2 n^2$', marker = '^')
    ax.plot(times, dOccTwoPh + 0.5, color='rosybrown', label='dOcc, $\omega_{\rm P}$ - $X^2 n$', marker = 'v')
    #ax.plot(times, NphOnePh, color='rosybrown', label='dOcc, $\omega_{\rm P} = 0.2$')

    plt.legend()

    # ax.set_ylim(-1e10, 1e10)

    plt.show()