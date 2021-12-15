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

    fileOnePh = h5py.File("../data/gsProp1PhGPH20NB20.hdf5", 'r')
    wPOnePh = fileOnePh['times'][()]
    dOccOnePh = fileOnePh['dOcc'][()]
    NptOnePh = fileOnePh['Npt'][()]
    NphOnePh = fileOnePh['N1ph'][()]

    fileTwoPh = h5py.File("../data/gsProp2PhGPH50NB4.hdf5", 'r')
    wPTwoPh = fileTwoPh['times'][()]
    dOccTwoPh = fileTwoPh['dOcc'][()]
    NptTwoPh = fileTwoPh['Npt'][()]
    NphTwoPh = fileTwoPh['N1ph'][()]

    fileTwoPhN6 = h5py.File("../data/gsProp2PhGPH50NB6.hdf5", 'r')
    dOccTwoPhN6 = fileTwoPh['dOcc'][()]

    fileTwoPhN8 = h5py.File("../data/gsProp2PhGPH50NB8.hdf5", 'r')
    dOccTwoPhN8 = fileTwoPh['dOcc'][()]

    fig = plt.figure()
    ax = fig.add_subplot(111)


    ax.plot(wPOnePh, dOccOnePh + 0.5, color='olive', label='dOcc, $\omega_{\rm P}$ - $X^2 n$', marker = '^')
    ax.plot(wPTwoPh, dOccTwoPh + 0.5, color='rosybrown', label='dOcc, $\omega_{\rm P}$ - $X^2 n^2$', marker = 'v')
    #ax.plot(wPTwoPh, NphTwoPh, color='olive', label=r'$\langle X^2 \rangle$ - $X^2 n$')
    #ax.plot(wPTwoPh, NptTwoPh, color='rosybrown', label=r'$\langle X^2 \rangle$ - $X^2 n^2$')



    plt.legend()

    # ax.set_ylim(-1e10, 1e10)

    plt.show()