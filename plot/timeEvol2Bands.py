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

    wD = .6

    # read in stuff
    file = h5py.File("./dataNB10/gsProp2BandsUa100Ub100Uud100Uss100epsA0epsB10gE20WD{}FD1000NB10TS20WPh5000.hdf5".format(int(np.rint(wD * 100))), 'r')

    #readInPrmsAndAssert(fileC, fileQ1)
    times = file['times'][()]
    pump = file['pump'][()]
    dOcc0 = file['dOcc0'][()]
    dOcc1 = file['dOcc1'][()]
    dOccUpDn = file['dOccUpDn'][()]
    dOccSigSig = file['dOccSigSig'][()]
    nPh = file['Nph1'][()]


    dOccIntra = dOcc0 + dOcc1
    dOccInter = dOccUpDn + dOccSigSig
    dOccTot = dOcc0 + dOcc1 + dOccUpDn + dOccSigSig

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3., 2.)

    #ax.plot(times, pump)
    ax.plot(times, dOccIntra, linewidth = 1., color = 'darkseagreen', label = r"$<n_{\alpha, \uparrow} n_{\alpha,\downarrow}>$")
    ax.plot(times, dOccInter, linewidth = 1., color = 'darkgreen', label = r"$<n_{\alpha, \sigma} n_{\beta,\sigma'}>$")
    #ax.plot(times, nPh, linewidth = 1., color = 'indianred', label = r"$N_{\rm ph}$")
    #ax.plot(times2, dOccTot, linewidth = .5)

    legend = ax.legend(fontsize=fontsize - 2, loc='lower left', bbox_to_anchor=(-0., 1.), edgecolor='black', ncol=2)
    legend.get_frame().set_alpha(1.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    plt.savefig('./savedPlots/timeEvol.png', format='png', bbox_inches='tight', dpi=600)



def plotTimeEvol2BandsManyPrms():
    print("plotting more beautiful time evolution")

    nWD = 60

    nPhMax = np.zeros(nWD)
    dOccIntraMax = np.zeros(nWD)
    dOccInterMax = np.zeros(nWD)

    wDArr = np.zeros(nWD)
    for wDInd, wD in enumerate(wDArr):
        wDArr[wDInd] = (wDInd + 1.) / 10.

    print(wDArr)

    for wDInd, wD in enumerate(wDArr):
        filename = "./dataNB10/gsProp2BandsUa100Ub100Uud100Uss100epsA0epsB10gE20WD{}FD1000NB10TS20WPh5000.hdf5".format(int(np.rint(wD * 100)))
        #print(filename)
        file = h5py.File(filename,'r')
        times = file['times'][()]
        pump = file['pump'][()]
        dOcc0 = file['dOcc0'][()]
        dOcc1 = file['dOcc1'][()]
        dOccUpDn = file['dOccUpDn'][()]
        dOccSigSig = file['dOccSigSig'][()]
        nPh = file['Nph1'][()]

        nPhMax[wDInd] = np.amax(nPh)
        dOccIntraMax[wDInd] = np.amin(dOcc0 + dOcc1)
        dOccInterMax[wDInd] = np.amax(dOccUpDn + dOccSigSig)

    #print(nPhMax)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3., 2.)

    ax.plot(wDArr, nPhMax, linewidth = 1., color = 'indianred', label = r'$N_{\rm ph}$')
    ax.plot(wDArr, 10 * dOccIntraMax, linewidth = 1., color = 'darkseagreen', label = r"$<n_{\alpha, \uparrow} n_{\alpha,\downarrow}>$")
    ax.plot(wDArr, 10 * dOccInterMax, linewidth = 1., color = 'darkgreen', label = r"$<n_{\alpha, \sigma} n_{\beta,\sigma'}>$")

    ax.axvline(1., linestyle = '--', linewidth = 0.8, color = 'gray')
    ax.axvline(5., linestyle = '--', linewidth = 0.8, color = 'gray')

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.3)

    ax.text(0.5, 2.5, r"$\varepsilon_1 {-} \varepsilon_0$", fontsize=10, bbox=boxProps)
    ax.text(4.7, 2.5, r"$\omega_{\rm ph}$", fontsize=10, bbox=boxProps)

    ax.set_ylim(0., 2.3)
    ax.set_xlim(0., 6.)

    ax.set_xlabel(r"$\omega_{\rm ph}$")
    ax.set_ylabel(r"$N$")

    legend = ax.legend(fontsize=fontsize, loc='upper left', bbox_to_anchor=(0., 1.), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(1.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    plt.savefig('./savedPlots/maxAsOfWD.png', format='png', bbox_inches='tight', dpi=600)


