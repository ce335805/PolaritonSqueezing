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

mpl.rcParams['font.family'] = 'Arial'
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

    wD = 9.

    # read in stuff
    file = h5py.File("../dataSmallTDriving/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB80gE40tH1WD{}FD1000NB6TS20WPh10000.hdf5".format(int(np.rint(wD * 100))), 'r')

    #readInPrmsAndAssert(fileC, fileQ1)
    times = file['times'][()]
    pump = file['pump'][()]
    dOcc0 = file['dOcc0'][()]
    dOcc1 = file['dOcc1'][()]
    dOccUpDn = file['dOccUpDn'][()]
    dOccSigSig = file['dOccSigSig'][()]
    n0 = file['n0'][()]
    n1 = file['n1'][()]
    nPh = file['Nph1'][()]


    dOccIntra = dOcc0 - 2. *  n0 * n0 + dOcc1 - 2. *  n1 * n1
    dOccInter = dOccUpDn - 2. * n0 * n1 + dOccSigSig - 2. *  n0 * n1

    #dOccIntra = n0
    #dOccInter = n1


    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3., 2.)

    #ax.plot(times, pump)
    ax.plot(times, dOccIntra, linewidth = 1., color = 'darkseagreen', label = r"$<n_{\alpha, \uparrow} n_{\alpha,\downarrow}>$")
    ax.plot(times, dOccInter, linewidth = 1., color = 'darkgreen', label = r"$<n_{\alpha, \sigma} n_{\beta,\sigma'}>$")
    #ax.plot(times, nPh, linewidth = 1., color = 'indianred', label = r"$N_{\rm ph}$")
    #ax.plot(times2, dOccTot, linewidth = .5)

    ax.set_yticks([0., 0.1, 0.2, 0.3, 0.4, 0.5])

    legend = ax.legend(fontsize=fontsize - 2, loc='lower left', bbox_to_anchor=(-0., 1.1), edgecolor='black', ncol=2)
    legend.get_frame().set_alpha(1.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    plt.savefig('./savedPlots/timeEvol.png', format='png', bbox_inches='tight', dpi=600)



def plotTimeEvol2BandsManyPrms():
    print("plotting more beautiful time evolution")

    nWD = 40

    nPhMax = np.zeros(nWD)
    dOccIntraMax = np.zeros(nWD)
    dOccInterMax = np.zeros(nWD)

    wDArr = np.zeros(nWD)
    for wDInd, wD in enumerate(wDArr):
        wDArr[wDInd] = (wDInd + 1.) / 4 + 6.

    print(wDArr)

    for wDInd, wD in enumerate(wDArr):
        #filename = "./dataNB10/gsProp2BandsUa100Ub100Uud100Uss100epsA0epsB10gE20WD{}FD1000NB10TS20WPh5000.hdf5".format(int(np.rint(wD * 100)))
        filename = "../dataTimeEvNb6/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB80gE40WD{}FD100NB6TS20WPh12000.hdf5".format(int(np.rint(wD * 100)))
        #print(filename)
        file = h5py.File(filename,'r')
        times = file['times'][()]
        pump = file['pump'][()]
        dOcc0 = file['dOcc0'][()]
        dOcc1 = file['dOcc1'][()]
        n0 = file['n0'][()]
        n1 = file['n1'][()]
        dOccUpDn = file['dOccUpDn'][()]
        dOccSigSig = file['dOccSigSig'][()]
        nPh = file['Nph1'][()]

        nPhMax[wDInd] = np.amax(nPh)
        dOccIntraMax[wDInd] = np.amax(dOcc0 - 2. * n0 * n0 + dOcc1 - 2. * n1 * n1)
        dOccInterMax[wDInd] = np.amin(dOccUpDn - 2. * n0 * n1 + dOccSigSig - 2. * n0 * n1)

    #print(nPhMax)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3., 2.)

    ax.plot(wDArr, nPhMax, linewidth = 1., color = 'indianred', label = r'$N_{\rm ph}$')
    #ax.plot(wDArr, dOccIntraMax, linewidth = 1., color = 'darkseagreen', label = r"$<n_{\alpha, \uparrow} n_{\alpha,\downarrow}>$")
    ax.plot(wDArr, dOccInterMax, linewidth = 1., color = 'darkgreen', label = r"$<n_{\alpha, \sigma} n_{\beta,\sigma'}>$")

    #ax.axvline(1., linestyle = '--', linewidth = 0.8, color = 'gray')
    #ax.axvline(5., linestyle = '--', linewidth = 0.8, color = 'gray')

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.3)

    #ax.text(0.5, 2.5, r"$\varepsilon_1 {-} \varepsilon_0$", fontsize=10, bbox=boxProps)
    #ax.text(4.7, 2.5, r"$\omega_{\rm ph}$", fontsize=10, bbox=boxProps)

    #ax.set_ylim(0., 2.3)
    #ax.set_xlim(0., 6.)

    ax.set_xlabel(r"$\omega_{\rm ph}$")
    ax.set_ylabel(r"$N$")

    #legend = ax.legend(fontsize=fontsize, loc='upper left', bbox_to_anchor=(0., 1.), edgecolor='black', ncol=1)
    #legend.get_frame().set_alpha(1.)
    #legend.get_frame().set_boxstyle('Square', pad=0.1)
    #legend.get_frame().set_linewidth(0.0)

    plt.savefig('./savedPlots/maxAsOfWD.png', format='png', bbox_inches='tight', dpi=600)


def plotTimeEvol2BandsManyWDWph():
    print("plotting more beautiful time evolution")


    wPhArr = np.array([8., 9., 9.8, 10.2, 11.5, 12.])
    #wPhArr = np.array([8., 9., 9.9, 10.1, 10.5, 11.5])
    wPh = len(wPhArr)
    nWD = 50


    dOccIntraMax = np.zeros((wPh, nWD))
    dOccInterMax = np.zeros((wPh, nWD))
    nPhMax = np.zeros((wPh, nWD))

    wDArr = np.zeros(nWD)
    for wDInd, wD in enumerate(wDArr):
        wDArr[wDInd] = (wDInd + 1.) / 5. + 5.

    for wPhInd, wPh in enumerate(wPhArr):
        for wDInd, wD in enumerate(wDArr):
            #filename = "./dataNB10/gsProp2BandsUa100Ub100Uud100Uss100epsA0epsB10gE20WD{}FD1000NB10TS20WPh5000.hdf5".format(int(np.rint(wD * 100)))
            filename = "../dataSmallTH/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE1tH1WD{}FD1000NB6TS20WPh{}.hdf5".format(int(np.rint(wD * 100)), int(np.rint(wPh * 1000)))
            #print(filename)
            file = h5py.File(filename,'r')
            times = file['times'][()]
            pump = file['pump'][()]
            dOcc0 = file['dOcc0'][()]
            dOcc1 = file['dOcc1'][()]
            n0 = file['n0'][()]
            n1 = file['n1'][()]
            dOccUpDn = file['dOccUpDn'][()]
            dOccSigSig = file['dOccSigSig'][()]
            nPh = file['Nph1'][()]

            nPhMax[wPhInd, wDInd] = np.amax(nPh)
            dOccIntraMax[wPhInd, wDInd] = np.amax(dOcc0 - 2. * n0 * n0 + dOcc1 - 2. * n1 * n1)
            dOccInterMax[wPhInd, wDInd] = np.amin(dOccUpDn - 2. * n0 * n1 + dOccSigSig - 2. * n0 * n1)

    #print(nPhMax)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3., 2.)

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')



    for wPhInd, wPh in enumerate(wPhArr):

        color = cmapPink(wPhInd / (len(wPhArr) + 1.))

        ax.plot(wDArr, nPhMax[wPhInd, :], linewidth = 1., color = color, label = r'$\omega = {}$'.format(wPh), zorder = 100)

        #ax.plot(wDArr, dOccIntraMax[wPhInd, :], linewidth = 1., color = color, label = r"$<n_{\alpha, \uparrow} n_{\alpha,\downarrow}>$", zorder = 100)
        #ax.plot(wDArr, dOccInterMax[wPhInd, :], linewidth = 1., color = color2, label = r"$<n_{\alpha, \sigma} n_{\beta,\sigma'}>$")

        ax.axvline(wPh, linestyle = '--', linewidth = 0.5, color = color)


    #ax.set_ylim(0., 0.135)
    ax.set_ylim(0., 1.6)
    ax.set_xlim(5.2, 15.)

    ax.set_xlabel(r"$\omega_{\rm Drive}$")
    ax.set_ylabel(r"$\rm{Max} \,\, N_{\rm plas}$")
    #ax.set_ylabel(r"$\rm{Max} \,\,$" + r"$\langle n_{\uparrow} n_{\downarrow} \rangle - \langle n \rangle^2$")

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.1)

    yPos = 0.14
    yPos = 1.65
    ax.text(7., yPos, r"$\omega_{\rm ph} = 8t_{\rm h}$", fontsize=8, bbox=boxProps, color = cmapPink(0. / (len(wPhArr) + 1.)))
    ax.text(9., yPos, r"$\omega_{\rm ph} = 9.8t_{\rm h}$", fontsize=8, bbox=boxProps, color = cmapPink(2. / (len(wPhArr) + 1.)))
    ax.text(11.5, yPos, r"$\omega_{\rm ph} = 12t_{\rm h}$", fontsize=8, bbox=boxProps, color = cmapPink(5. / (len(wPhArr) + 1.)))
    #ax.text(4.7, 2.5, r"$\omega_{\rm ph}$", fontsize=10, bbox=boxProps)


    #legend = ax.legend(fontsize=fontsize, loc='upper left', bbox_to_anchor=(0., 1.), edgecolor='black', ncol=1)
    #legend.get_frame().set_alpha(1.)
    #legend.get_frame().set_boxstyle('Square', pad=0.1)
    #legend.get_frame().set_linewidth(0.0)

    plt.savefig('./savedPlots/maxNAsOfWDAndWPh.png', format='png', bbox_inches='tight', dpi=600)
