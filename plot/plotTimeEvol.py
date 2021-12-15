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

def avOverArr(arr, avOver):

    averagedArr = np.zeros(arr.shape)

    for ind in np.arange(len(arr)):
        for avInd in np.arange(avOver):
            averagedArr[ind] += arr[(ind - (avInd - avOver // 2)) % len(arr)] / avOver
    return averagedArr




def plotQCCompareNSqr():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileC = h5py.File("../data/tEvol1PhGPH20WP0WD200FD200NB10.hdf5", 'r')
    fileQ1 = h5py.File("../data/tEvol1PhGPH20WP50WD227FD200NB10.hdf5", 'r')
    fileQ2 = h5py.File("../data/tEvol1PhGPH20WP50WD177FD200NB10.hdf5", 'r')

    #readInPrmsAndAssert(fileC, fileQ1)
    wPhQ = (fileQ1['wPh'][()])[0]
    wPtQ = (fileQ1['wPt'][()])[0]
    wPQ = (fileQ1['wP'][()])[0]

    wDriveC = (fileC['wDrive'][()])[0]
    wDriveQ = (fileQ1['wDrive'][()])[0]


    wMinus = plotPolaritonFreqs.calcWMinus(wPhQ, wPtQ, wPQ)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhQ, wPtQ, wPQ)

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDriveQ))


    times = fileC['times'][()]
    times = times / (2. * np.pi / wDriveC)
    pump = fileC['pump'][()]

    dOccC = fileC['dOcc'][()]
    dOccQ = fileQ1['dOcc'][()]
    dOccQ2 = fileQ2['dOcc'][()]

    Xpt = fileC['Xpt'][()]
    XptSqr = fileC['XptSqr'][()]
    Npt = fileC['Npt'][()]
    Xph = fileC['X1ph'][()]
    XphSqr = fileC['X1phSqr'][()]
    Nph = fileC['N1ph'][()]

    NptQ = fileQ1['Npt'][()]
    NphQ = fileQ1['N1ph'][()]
    XphSqrQ = fileQ1['X1phSqr'][()]

    NptQ2 = fileQ2['Npt'][()]
    NphQ2 = fileQ2['N1ph'][()]
    XphSqrQ2 = fileQ2['X1phSqr'][()]


    fig = plt.figure()
    fig.set_size_inches(3., 3.)

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex = ax1)
    ax3 = fig.add_subplot(gs[2], sharex = ax1)
    ax4 = fig.add_subplot(gs[3], sharex = ax1)

    ax1.tick_params(direction='inout', length=6, width=.5)
    ax2.tick_params(direction='inout', length=6, width=.5)
    ax3.tick_params(direction='inout', length=6, width=.5)

    #ax1.plot(times, Nph * 5., color='rosybrown', label=r'$N_{ph} \times 5$', linewidth = 1.5)
    #ax1.plot(times, XphSqr * 0.2, color='olive', label=r'$\langle X^2 \rangle \times \omega_{\rm ph}$', linewidth = 1.)
    ax1.plot(times, dOccC + 0.5, color='olive', label='Phonon Driving', linewidth = 1.)
    ax1.plot(times, dOccQ + 0.5, color='rosybrown', label='Up-Polariton Driving', linewidth = 1.)
    ax1.plot(times, dOccQ2 + 0.5, color='gray', label='Low-Polariton Driving', linewidth = 1.)

    ax2.plot(times, Nph, color='olive', label='Classical Driving', linewidth = 1.)
    ax2.plot(times, NphQ, color='rosybrown', label='Driven Cavity', linewidth = 1.)
    ax2.plot(times, NphQ2, color='gray', label='Driven Cavity', linewidth = 1.)

    ax3.plot(times, Npt, color='olive', label='Classical Driving', linewidth = 1.)
    ax3.plot(times, NptQ, color='rosybrown', label='Driven Cavity', linewidth = 1.)
    ax3.plot(times, NptQ2, color='gray', label='Driven Cavity', linewidth = 1.)

    ax4.plot(times, pump, color='cornflowerblue', label='pump', linewidth = 1.)

    #ax3.plot(times, XphSqr, color='olive', label='Driven Cavity', linewidth = 1.)
    #ax3.plot(times, XphSqrQ, color='rosybrown', label='Driven Cavity', linewidth = 1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize = fontsize)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize = fontsize)
    ax4.set_ylabel(r"$F(t)$", fontsize = fontsize)


    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize = fontsize)


    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim(3., 13)
    ax2.set_xlim(3., 13)
    ax3.set_xlim(3., 13)
    ax4.set_xlim(3., 13)
    ax1.set_xticks([4, 6, 8, 10, 12])
    ax2.set_xticks([4, 6, 8, 10, 12])
    ax3.set_xticks([4, 6, 8, 10, 12])
    ax3.set_xticklabels(["$4$", "$6$", "$8$", "$10$", "$12$"])


    ax1.set_ylim(0.08, 0.21)

    ax1.set_yticks([0.1, 0.2])
    ax1.set_yticklabels(["$0.1$", "$0.2$"])

    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(["$0$", "$1$"])

    ax3.set_yticks([0, 1])
    ax3.set_yticklabels(["$0$", "$1$"])

    ax4.set_yticks([-0.2, 0., 0.2])
    ax4.set_yticklabels(["$-0.2$", "$0$", "$0.2$"])


    legend1 = ax1.legend(fontsize = fontsize - 4, loc = 'upper left', bbox_to_anchor=(.0, 1.), edgecolor = 'black', ncol = 1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(.0, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    #plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    #plt.show()

    plt.savefig('testPlot.png', format='png', bbox_inches='tight', dpi = 600)


def plotAsOfWP():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../data/tEvol1PhGPH20WP0WD200FD200NB10.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol1PhGPH20WP10WD200FD200NB10.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol1PhGPH20WP20WD200FD200NB10.hdf5", 'r')
    fileWP4 = h5py.File("../data/tEvol1PhGPH20WP50WD200FD200NB10.hdf5", 'r')

    #readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP4 = (fileWP4['wPh'][()])[0]
    wPtWP4 = (fileWP4['wPt'][()])[0]
    wPWP4 = (fileWP4['wP'][()])[0]




    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP4, wPtWP4, wPWP4)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP4, wPtWP4, wPWP4)

    wDrive = (fileWP1['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDrive))


    times = fileWP1['times'][()]
    times = times / (2. * np.pi / wDrive)
    pump = fileWP1['pump'][()]

    dOccWP1 = fileWP1['dOcc'][()]
    dOccWP2 = fileWP2['dOcc'][()]
    dOccWP3 = fileWP3['dOcc'][()]
    dOccWP200 = fileWP4['dOcc'][()]

    NptWP1 = fileWP1['Npt'][()]
    NptWP2 = fileWP2['Npt'][()]
    NptWP3 = fileWP3['Npt'][()]
    NptWP4 = fileWP4['Npt'][()]



    NphWP1 = fileWP1['N1ph'][()]
    NphWP2 = fileWP2['N1ph'][()]
    NphWP3 = fileWP3['N1ph'][()]
    NphWP4 = fileWP4['N1ph'][()]


    fig = plt.figure()
    fig.set_size_inches(3., 3.)

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex = ax1)
    ax3 = fig.add_subplot(gs[2], sharex = ax1)
    ax4 = fig.add_subplot(gs[3], sharex = ax1)

    ax1.tick_params(direction='inout', length=6, width=.5)
    ax2.tick_params(direction='inout', length=6, width=.5)
    ax3.tick_params(direction='inout', length=6, width=.5)

    linewidth = 0.8

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#FF785A'
    colorPump = '#5295A4'

    ax1.plot(times, dOccWP1 + 0.5, color=color1, label=r'$\omega_{\rm P} = 0$', linewidth = linewidth)
    ax1.plot(times, dOccWP2 + 0.5, color=color2, label=r'$\omega_{\rm P} = 0.1$', linewidth = linewidth)
    ax1.plot(times, dOccWP3 + 0.5, color=color3, label=r'$\omega_{\rm P} = 0.2$', linewidth = linewidth)
    ax1.plot(times, dOccWP200 + 0.5, color=color4, label=r'$\omega_{\rm P} = 0.5$', linewidth = linewidth)

    ax2.plot(times, NphWP1, color=color1, linewidth = linewidth)
    ax2.plot(times, NphWP2, color=color2, linewidth = linewidth)
    ax2.plot(times, NphWP3, color=color3, linewidth = linewidth)
    ax2.plot(times, NphWP4, color=color4, linewidth = linewidth)

    ax3.plot(times, NptWP1, color=color1, linewidth = linewidth)
    ax3.plot(times, NptWP2, color=color2, linewidth = linewidth)
    ax3.plot(times, NptWP3, color=color3, linewidth = linewidth)
    ax3.plot(times, NptWP4, color=color4, linewidth = linewidth)

    ax4.plot(times, pump, color=colorPump, label='pump', linewidth = 1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize = fontsize)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize = fontsize)
    ax4.set_ylabel(r"$F(t)$", fontsize = fontsize)


    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize = fontsize)


    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim(3., 17)
    ax2.set_xlim(3., 17)
    ax3.set_xlim(3., 17)
    ax4.set_xlim(3., 17)
    ax1.set_xticks([4, 6, 8, 10, 12])
    ax2.set_xticks([4, 6, 8, 10, 12])
    ax3.set_xticks([4, 6, 8, 10, 12])
    ax3.set_xticklabels(["$4$", "$6$", "$8$", "$10$", "$12$"])


    #ax1.set_ylim(0.08, 0.21)
    ax1.set_yticks([0.1, 0.13, 0.16])
    ax1.set_yticklabels(["$0.1$", "$0.13$", "$0.16$"])

    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(["$0$", "$1$"])

    ax3.set_yticks([0, 1, 2])
    ax3.set_yticklabels(["$0$", "$1$", "$2$"])

    ax4.set_yticks([-0.2, 0., 0.2])
    ax4.set_yticklabels(["$-0.2$", "$0$", "$0.2$"])


    legend1 = ax1.legend(fontsize = fontsize - 4, loc = 'upper left', bbox_to_anchor=(.0, 1.), edgecolor = 'black', ncol = 1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(0.7, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    #plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    #plt.show()

    plt.savefig('functionOfWP.png', format='png', bbox_inches='tight', dpi = 600)


def plotAsOfWPLin():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../data/tEvol2PhGPH50WP0WD280FD200NB4.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol2PhGPH50WP10WD280FD200NB4.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol2PhGPH50WP20WD280FD200NB4.hdf5", 'r')
    fileWP4 = h5py.File("../data/tEvol2PhGPH50WP50WD280FD200NB4.hdf5", 'r')

    #readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP4 = (fileWP4['wPh'][()])[0]
    wPtWP4 = (fileWP4['wPt'][()])[0]
    wPWP4 = (fileWP4['wP'][()])[0]




    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP4, wPtWP4, wPWP4)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP4, wPtWP4, wPWP4)

    wDrive = (fileWP1['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDrive))


    times = fileWP1['times'][()]
    times = times / (2. * np.pi / wDrive)
    pump = fileWP1['pump'][()]

    dOccWP1 = fileWP1['dOcc'][()]
    dOccWP2 = fileWP2['dOcc'][()]
    dOccWP3 = fileWP3['dOcc'][()]
    dOccWP200 = fileWP4['dOcc'][()]

    NptWP1 = fileWP1['Npt'][()]
    NptWP2 = fileWP2['Npt'][()]
    NptWP3 = fileWP3['Npt'][()]
    NptWP4 = fileWP4['Npt'][()]



    NphWP1 = fileWP1['N1ph'][()]
    NphWP2 = fileWP2['N1ph'][()]
    NphWP3 = fileWP3['N1ph'][()]
    NphWP4 = fileWP4['N1ph'][()]


    fig = plt.figure()
    fig.set_size_inches(3., 3.)

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex = ax1)
    ax3 = fig.add_subplot(gs[2], sharex = ax1)
    ax4 = fig.add_subplot(gs[3], sharex = ax1)

    ax1.tick_params(direction='inout', length=6, width=.5)
    ax2.tick_params(direction='inout', length=6, width=.5)
    ax3.tick_params(direction='inout', length=6, width=.5)

    linewidth = .8

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#FF785A'
    colorPump = '#5295A4'

    ax1.plot(times, dOccWP1 + 0.5, color=color1, label=r'$\omega_{\rm P} = 0$', linewidth = linewidth, zorder = 4)
    ax1.plot(times, dOccWP2 + 0.5, color=color2, label=r'$\omega_{\rm P} = 0.1$', linewidth = linewidth, zorder = 3)
    ax1.plot(times, dOccWP3 + 0.5, color=color3, label=r'$\omega_{\rm P} = 0.2$', linewidth = linewidth, zorder = 2)
    ax1.plot(times, dOccWP200 + 0.5, color=color4, label=r'$\omega_{\rm P} = 0.5$', linewidth = linewidth, zorder = 1)

    ax2.plot(times, NphWP1, color=color1, linewidth = linewidth)
    ax2.plot(times, NphWP2, color=color2, linewidth = linewidth)
    ax2.plot(times, NphWP3, color=color3, linewidth = linewidth)
    ax2.plot(times, NphWP4, color=color4, linewidth = linewidth)

    ax3.plot(times, NptWP1, color=color1, linewidth = linewidth)
    ax3.plot(times, NptWP2, color=color2, linewidth = linewidth)
    ax3.plot(times, NptWP3, color=color3, linewidth = linewidth)
    ax3.plot(times, NptWP4, color=color4, linewidth = linewidth)

    ax4.plot(times, pump, color=colorPump, label='pump', linewidth = 1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize = fontsize)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize = fontsize)
    ax4.set_ylabel(r"$F(t)$", fontsize = fontsize)


    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize = fontsize)


    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim(3., 17)
    ax2.set_xlim(3., 17)
    ax3.set_xlim(3., 17)
    ax4.set_xlim(3., 17)
    ax1.set_xticks([4, 6, 8, 10, 12])
    ax2.set_xticks([4, 6, 8, 10, 12])
    ax3.set_xticks([4, 6, 8, 10, 12])
    ax3.set_xticklabels(["$4$", "$6$", "$8$", "$10$", "$12$"])


    #ax1.set_ylim(0.08, 0.21)
    #ax1.set_yticks([0.1, 0.13, 0.16])
    #ax1.set_yticklabels(["$0.1$", "$0.13$", "$0.16$"])
#
    #ax2.set_yticks([0, 1])
    #ax2.set_yticklabels(["$0$", "$1$"])
#
    #ax3.set_yticks([0, 1, 2])
    #ax3.set_yticklabels(["$0$", "$1$", "$2$"])
#
    #ax4.set_yticks([-0.2, 0., 0.2])
    #ax4.set_yticklabels(["$-0.2$", "$0$", "$0.2$"])


    legend1 = ax1.legend(fontsize = fontsize - 4, loc = 'upper left', bbox_to_anchor=(.0, 1.), edgecolor = 'black', ncol = 1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(0.7, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    #plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    #plt.show()

    plt.savefig('functionOfWPLin.png', format='png', bbox_inches='tight', dpi = 600)


def readInPrmsAndAssert(fileC, fileQ):
    wPhC = (fileC['wPh'][()])[0]
    wPhQ = (fileQ['wPh'][()])[0]
    assert (wPhC == wPhQ)
    wPtC = (fileC['wPt'][()])[0]
    wPtQ = (fileQ['wPt'][()])[0]
    assert (wPtC == wPtQ)
    tHopC = (fileC['tHop'][()])[0]
    tHopQ = (fileQ['tHop'][()])[0]
    assert (tHopQ == tHopC)
    UC = (fileC['U'][()])[0]
    UQ = (fileQ['U'][()])[0]
    assert (UC == UQ)
    gPhC = (fileC['gPh'][()])[0]
    gPhQ = (fileQ['gPh'][()])[0]
    assert (gPhQ == gPhC)
    wPC = (fileC['wP'][()])[0]
    wPQ = (fileQ['wP'][()])[0]
    assert (wPC == 0.)
    assert (wPQ != 0.)

    wDriveC = (fileC['wDrive'][()])[0]
    wDriveQ = (fileQ['wDrive'][()])[0]
    #assert (wDriveQ == wDriveC)
    fDriveC = (fileC['fDrive'][()])[0]
    fDriveQ = (fileQ['fDrive'][()])[0]
    assert (fDriveQ == fDriveC)
    tsPerDrivePeriodC = (fileC['timePointsPerDrivingPeriod'][()])[0]
    tsPerDrivePeriodQ = (fileQ['timePointsPerDrivingPeriod'][()])[0]
    assert (tsPerDrivePeriodQ == tsPerDrivePeriodC)

    nPhononC = (fileC['dimPhonon'][()])[0]
    nPhononQ = (fileQ['dimPhonon'][()])[0]
    assert (nPhononQ == nPhononC)
    nPhotonC = (fileC['dimPhoton'][()])[0]
    assert (nPhotonC == 1)
    nPhotonQ = (fileQ['dimPhoton'][()])[0]
    assert (nPhotonQ == nPhononQ)

    print("tHop = {}".format(tHopC))
    print("U = {}".format(UC))
    print("wPh = {}".format(wPhC))
    print("wPt = {}".format(wPtC))
    print("gPh = {}".format(gPhC))
    print("wP = {}".format(wPC))
    print("wDrive = {}".format(wDriveC))
    print("fDrive = {}".format(fDriveC))
    print("time points per driving period = {}".format(tsPerDrivePeriodC))
    print("n-Phonon = {}".format(nPhononC))
    print("n-Photon = {}".format(nPhotonC))