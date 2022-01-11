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
    fileWP1 = h5py.File("../data/tEvol2PhGPH20WP71WD200FD200NB6TS80.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol2PhGPH20WP71WD156FD200NB6TS80.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol2PhGPH20WP71WD256FD200NB6TS80.hdf5", 'r')
    fileWP4 = h5py.File("../data/tEvol2PhGPH20WP71WD256FD200NB6TS80.hdf5", 'r')

    #fileWP1 = h5py.File("../data/tEvol2PhGPH20WP0WD200FD200NB4TS80.hdf5", 'r')
    #fileWP2 = h5py.File("../data/tEvol2PhGPH20WP7WD200FD200NB10TS80.hdf5", 'r')
    #fileWP3 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD200NB10TS80.hdf5", 'r')
    #fileWP4 = h5py.File("../data/tEvol2PhGPH20WP28WD200FD200NB10TS80.hdf5", 'r')

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

    times = fileWP2['times'][()]
    times = times / (2. * np.pi / wDrive)
    pump = fileWP2['pump'][()]

    dOccWP1 = fileWP1['dOcc'][()]
    dOccWP2 = fileWP2['dOcc'][()]
    dOccWP3 = fileWP3['dOcc'][()]
    dOccWP4 = fileWP4['dOcc'][()]

    NptWP1 = fileWP1['Npt'][()]
    NptWP2 = fileWP2['Npt'][()]
    NptWP3 = fileWP3['Npt'][()]
    NptWP4 = fileWP4['Npt'][()]

    NphWP1 = fileWP1['N1ph'][()]
    NphWP2 = fileWP2['N1ph'][()]
    NphWP3 = fileWP3['N1ph'][()]
    NphWP4 = fileWP4['N1ph'][()]

    fig = plt.figure()
    fig.set_size_inches(6., 3.)

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
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax1.plot(times, dOccWP1 + 0.5, color=color1, label=r'$\omega_{\rm P} = 0$', linewidth = linewidth)
    ax1.plot(times, dOccWP2 + 0.5, color=color2, label=r'$\omega_{\rm P} = 0.1$', linewidth = linewidth)
    ax1.plot(times, dOccWP3 + 0.5, color=color3, label=r'$\omega_{\rm P} = 0.2$', linewidth = linewidth)
    ax1.plot(times, dOccWP4 + 0.5, color=color4, label=r'$\omega_{\rm P} = 0.5$', linewidth = linewidth)

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

    ax1.set_xlim(3., 19)
    ax2.set_xlim(3., 19)
    ax3.set_xlim(3., 19)
    ax4.set_xlim(3., 19)
    ax1.set_xticks([4, 8, 12, 16])
    ax2.set_xticks([4, 8, 12, 16])
    ax3.set_xticks([4, 8, 12, 16])
    ax3.set_xticklabels(["$4$", "$8$", "$12$", "$16$"])


    #ax1.set_ylim(0.08, 0.21)
    #ax1.set_yticks([0.1, 0.15, 0.2])
    #ax1.set_yticklabels(["$0.1$", "$0.15$", "$0.2$"])
#
    ax2.set_yticks([0, 0.5, 1])
    ax2.set_yticklabels(["$0$", "$0.5$", '$1$'])
#
    ax3.set_yticks([0, 1, 2])
    ax3.set_yticklabels(["$0$", "$1$", "$2$"])
#
    ax4.set_ylim(-0.4, 0.4)
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

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    #plt.savefig('functionOfWP.png', format='png', bbox_inches='tight', dpi = 600)


def plotAsOfWPUSC():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../data/tEvol2PhGPH20WP71WD200FD200NB6TS80.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol2PhGPH20WP71WD156FD200NB6TS80.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol2PhGPH20WP71WD256FD200NB6TS80.hdf5", 'r')

    #readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP1 = (fileWP1['wPh'][()])[0]
    wPtWP1 = (fileWP1['wPt'][()])[0]
    wPWP1 = (fileWP1['wP'][()])[0]
    
    wPhWP2 = (fileWP2['wPh'][()])[0]
    wPtWP2 = (fileWP2['wPt'][()])[0]
    wPWP2 = (fileWP2['wP'][()])[0]
    
    wPhWP3 = (fileWP3['wPh'][()])[0]
    wPtWP3 = (fileWP3['wPt'][()])[0]
    wPWP3 = (fileWP3['wP'][()])[0]

    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP1, wPtWP1, wPWP1)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP1, wPtWP1, wPWP1)

    wDriveWP1 = (fileWP1['wDrive'][()])[0]
    wDriveWP2 = (fileWP2['wDrive'][()])[0]
    wDriveWP3 = (fileWP3['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDriveWP1))

    timesWP1 = fileWP1['times'][()]
    #timesWP1 = timesWP1 / (2. * np.pi / wDriveWP1)
    timesWP1 = timesWP1 / (2. * np.pi / 2.)
    pumpWP1 = fileWP1['pump'][()]

    timesWP2 = fileWP2['times'][()]
    #timesWP2 = timesWP2 / (2. * np.pi / wDriveWP2)
    timesWP2 = timesWP2 / (2. * np.pi / 2.)
    pumpWP2 = fileWP2['pump'][()]

    timesWP3 = fileWP3['times'][()]
    #timesWP3 = timesWP3 / (2. * np.pi / wDriveWP3)
    timesWP3 = timesWP3 / (2. * np.pi / 2.)
    pumpWP3 = fileWP3['pump'][()]

    dOccWP1 = fileWP1['dOcc'][()]
    dOccWP2 = fileWP2['dOcc'][()]
    dOccWP3 = fileWP3['dOcc'][()]

    NptWP1 = fileWP1['Npt'][()]
    NptWP2 = fileWP2['Npt'][()]
    NptWP3 = fileWP3['Npt'][()]

    NphWP1 = fileWP1['N1ph'][()]
    NphWP2 = fileWP2['N1ph'][()]
    NphWP3 = fileWP3['N1ph'][()]

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
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax1.plot(timesWP1, dOccWP1 + 0.5, color=color1, label=r'$\omega_{\rm P} = 0$', linewidth = linewidth)
    ax1.plot(timesWP2, dOccWP2 + 0.5, color=color2, label=r'$\omega_{\rm P} = 0.1$', linewidth = linewidth)
    ax1.plot(timesWP3, dOccWP3 + 0.5, color=color3, label=r'$\omega_{\rm P} = 0.2$', linewidth = linewidth)

    ax2.plot(timesWP1, NphWP1, color=color1, linewidth = linewidth)
    ax2.plot(timesWP2, NphWP2, color=color2, linewidth = linewidth)
    ax2.plot(timesWP3, NphWP3, color=color3, linewidth = linewidth)

    ax3.plot(timesWP1, NptWP1, color=color1, linewidth = linewidth)
    ax3.plot(timesWP2, NptWP2, color=color2, linewidth = linewidth)
    ax3.plot(timesWP3, NptWP3, color=color3, linewidth = linewidth)

    ax4.plot(timesWP1, pumpWP1, color=color1, label='pump', linewidth = 1.)
    ax4.plot(timesWP2, pumpWP2, color=color2, label='pump', linewidth = 1.)
    ax4.plot(timesWP3, pumpWP3, color=color3, label='pump', linewidth = 1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize = fontsize)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize = fontsize)
    ax4.set_ylabel(r"$F(t)$", fontsize = fontsize)


    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize = fontsize)


    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)

    ax1.set_xlim(1., 15)
    ax2.set_xlim(1., 15)
    ax3.set_xlim(1., 15)
    ax4.set_xlim(1., 15)
    ax4.set_xticks([4, 8, 12])
    ax4.set_xticklabels(["$4$", "$8$", "$12$"])


    #ax1.set_ylim(0.08, 0.21)
    #ax1.set_yticks([0.1, 0.15, 0.2])
    #ax1.set_yticklabels(["$0.1$", "$0.15$", "$0.2$"])
#
    ax2.set_yticks([0, 0.5, 1])
    ax2.set_yticklabels(["$0$", "$0.5$", '$1$'])
#
    ax3.set_yticks([0, 1, 2])
    ax3.set_yticklabels(["$0$", "$1$", "$2$"])
#
    ax4.set_ylim(-0.4, 0.4)
    ax4.set_yticks([-0.2, 0., 0.2])
    ax4.set_yticklabels(["$-0.2$", "$0$", "$0.2$"])


    legend1 = ax1.legend(fontsize = fontsize - 4, loc = 'upper left', bbox_to_anchor=(.0, 1.), edgecolor = 'black', ncol = 1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    #legend4 = ax4.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(0.7, 1.1), edgecolor='black', ncol=1)
    #legend4.get_frame().set_alpha(0.)
    #legend4.get_frame().set_boxstyle('Square', pad=0.1)
    #legend4.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    #plt.savefig('functionOfWP.png', format='png', bbox_inches='tight', dpi = 600)


def plotQuadUSCpSC():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../data/tEvol2PhGPH20WP71WD200FD300NB10TS80WPh2022.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol2PhGPH20WP71WD156FD300NB10TS80WPh2022.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol2PhGPH20WP71WD256FD300NB10TS80WPh2022.hdf5", 'r')

    fileWP4 = h5py.File("../data/tEvol2PhGPH20WP4WD200FD300NB10TS80WPh2022.hdf5", 'r')
    fileWP5 = h5py.File("../data/tEvol2PhGPH20WP7WD200FD300NB10TS80WPh2022.hdf5", 'r')

    # readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP1 = (fileWP1['wPh'][()])[0]
    wPtWP1 = (fileWP1['wPt'][()])[0]
    wPWP1 = (fileWP1['wP'][()])[0]

    wPhWP2 = (fileWP2['wPh'][()])[0]
    wPtWP2 = (fileWP2['wPt'][()])[0]
    wPWP2 = (fileWP2['wP'][()])[0]

    wPhWP3 = (fileWP3['wPh'][()])[0]
    wPtWP3 = (fileWP3['wPt'][()])[0]
    wPWP3 = (fileWP3['wP'][()])[0]
    
    wPhWP4 = (fileWP4['wPh'][()])[0]
    wPtWP4 = (fileWP4['wPt'][()])[0]
    wPWP4 = (fileWP4['wP'][()])[0]
    
    wPhWP5 = (fileWP5['wPh'][()])[0]
    wPtWP5 = (fileWP5['wPt'][()])[0]
    wPWP5 = (fileWP5['wP'][()])[0]

    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP1, wPtWP1, wPWP1)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP1, wPtWP1, wPWP1)

    wDriveWP1 = (fileWP1['wDrive'][()])[0]
    wDriveWP2 = (fileWP2['wDrive'][()])[0]
    wDriveWP3 = (fileWP3['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDriveWP1))

    timesWP1 = fileWP1['times'][()]
    # timesWP1 = timesWP1 / (2. * np.pi / wDriveWP1)
    timesWP1 = timesWP1 / (2. * np.pi / 2.)
    pumpWP1 = fileWP1['pump'][()]

    timesWP2 = fileWP2['times'][()]
    # timesWP2 = timesWP2 / (2. * np.pi / wDriveWP2)
    timesWP2 = timesWP2 / (2. * np.pi / 2.)
    pumpWP2 = fileWP2['pump'][()]

    timesWP3 = fileWP3['times'][()]
    # timesWP3 = timesWP3 / (2. * np.pi / wDriveWP3)
    timesWP3 = timesWP3 / (2. * np.pi / 2.)
    pumpWP3 = fileWP3['pump'][()]
    
    timesWP4 = fileWP4['times'][()]
    timesWP4 = timesWP4 / (2. * np.pi / 2.)
    pumpWP4 = fileWP4['pump'][()]
    
    timesWP5 = fileWP5['times'][()]
    timesWP5 = timesWP5 / (2. * np.pi / 2.)
    pumpWP5 = fileWP5['pump'][()]

    dOccWP1 = fileWP1['dOcc'][()]
    dOccWP2 = fileWP2['dOcc'][()]
    dOccWP3 = fileWP3['dOcc'][()]
    dOccWP4 = fileWP4['dOcc'][()]
    dOccWP5 = fileWP5['dOcc'][()]

    NptWP1 = fileWP1['Npt'][()]
    NptWP2 = fileWP2['Npt'][()]
    NptWP3 = fileWP3['Npt'][()]
    NptWP4 = fileWP4['Npt'][()]
    NptWP5 = fileWP5['Npt'][()]

    NphWP1 = fileWP1['N1ph'][()]
    NphWP2 = fileWP2['N1ph'][()]
    NphWP3 = fileWP3['N1ph'][()]
    NphWP4 = fileWP4['N1ph'][()]
    NphWP5 = fileWP5['N1ph'][()]

    fig = plt.figure()
    fig.set_size_inches(6., 3.)

    gs = gridspec.GridSpec(4, 2, height_ratios=[2, 1, 1, 1])
    ax11 = fig.add_subplot(gs[0, 1])
    ax21 = fig.add_subplot(gs[1, 1], sharex=ax11)
    ax31 = fig.add_subplot(gs[2, 1], sharex=ax11)
    ax41 = fig.add_subplot(gs[3, 1], sharex=ax11)

    ax12 = fig.add_subplot(gs[0, 0], sharey = ax11)
    ax22 = fig.add_subplot(gs[1, 0], sharex=ax12, sharey = ax21)
    ax32 = fig.add_subplot(gs[2, 0], sharex=ax12, sharey = ax31)
    ax42 = fig.add_subplot(gs[3, 0], sharex=ax12, sharey = ax41)

    ax11.tick_params(direction='inout', length=4, width=.8)
    ax21.tick_params(direction='inout', length=4, width=.8)
    ax31.tick_params(direction='inout', length=4, width=.8)
    ax41.tick_params(direction='inout', length=4, width=.8)

    ax12.tick_params(direction='inout', length=4, width=.8)
    ax22.tick_params(direction='inout', length=4, width=.8)
    ax32.tick_params(direction='inout', length=4, width=.8)
    ax42.tick_params(direction='inout', length=4, width=.8)

    ax11.yaxis.tick_right()
    ax21.yaxis.tick_right()
    ax31.yaxis.tick_right()
    ax41.yaxis.tick_right()

    linewidth = 0.8

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax11.plot(timesWP1, dOccWP1 + 0.5, color=color1, label=r'$\omega_{\rm Drive} = \omega_{\rm phot}$', linewidth=linewidth, zorder = 10)
    ax11.plot(timesWP2, dOccWP2 + 0.5, color=color2, label=r'$\omega_{\rm Drive} = \omega_-$', linewidth=linewidth, zorder = 5)
    ax11.plot(timesWP3, dOccWP3 + 0.5, color=color3, label=r'$\omega_{\rm Drive} = \omega_+$', linewidth=linewidth)

    ax21.plot(timesWP1, 2. * NphWP1, color=color1, linewidth=linewidth, zorder = 10)
    ax21.plot(timesWP2, 2. * NphWP2, color=color2, linewidth=linewidth, zorder = 5)
    ax21.plot(timesWP3, 2. * NphWP3, color=color3, linewidth=linewidth)

    ax31.plot(timesWP1, NptWP1, color=color1, linewidth=linewidth, zorder = 10)
    ax31.plot(timesWP2, NptWP2, color=color2, linewidth=linewidth, zorder = 5)
    ax31.plot(timesWP3, NptWP3, color=color3, linewidth=linewidth)

    ax41.plot(timesWP1, pumpWP1 * 2. / 3., color=color1, label='pump', linewidth=1., zorder = 10)
    ax41.plot(timesWP2, pumpWP2 * 2. / 3., color=color2, label='pump', linewidth=1., zorder = 5)
    ax41.plot(timesWP3, pumpWP3 * 2. / 3., color=color3, label='pump', linewidth=1.)

    ax12.plot(timesWP4, dOccWP4 + 0.5, color=color1, label=r'$\omega_{\rm phot} = 20 \omega_{\rm P}$', linewidth=linewidth)
    ax12.plot(timesWP5, dOccWP5 + 0.5, color=color3, label=r'$\omega_{\rm phot} = 40 \omega_{\rm P}$', linewidth=linewidth)

    ax22.plot(timesWP4, 2. * NphWP4, color=color1, linewidth=linewidth)
    ax22.plot(timesWP5, 2. * NphWP5, color=color3, linewidth=linewidth)

    ax32.plot(timesWP4, NptWP4, color=color1, linewidth=linewidth)
    ax32.plot(timesWP5, NptWP5, color=color3, linewidth=linewidth)

    ax42.plot(timesWP4, pumpWP4 * 2. / 3., color=color1, label='pump', linewidth=1.)
    ax42.plot(timesWP5, pumpWP5 * 2. / 3., color=color3, label='pump', linewidth=1., linestyle = '--', dashes = [4, 4])

    ax12.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize, labelpad=3.)
    ax22.set_ylabel(r"$N_{\rm phon}$", fontsize=fontsize, labelpad=15.)
    ax32.set_ylabel(r"$N_{\rm phot}$", fontsize=fontsize, labelpad=15.)
    ax42.set_ylabel(r"$F(t) / F_0$", fontsize=fontsize, labelpad=0.)

    ax41.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm phot}]$", fontsize=fontsize)
    ax42.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm phot}]$", fontsize=fontsize)

    plt.setp(ax11.get_xticklabels(), visible=False)
    plt.setp(ax21.get_xticklabels(), visible=False)
    plt.setp(ax31.get_xticklabels(), visible=False)

    plt.setp(ax12.get_xticklabels(), visible=False)
    plt.setp(ax22.get_xticklabels(), visible=False)
    plt.setp(ax32.get_xticklabels(), visible=False)

    plt.setp(ax11.get_yticklabels(), visible=False)
    plt.setp(ax21.get_yticklabels(), visible=False)
    plt.setp(ax31.get_yticklabels(), visible=False)
    plt.setp(ax41.get_yticklabels(), visible=False)


    ax11.set_xlim(3., 18)
    ax21.set_xlim(3., 18)
    ax31.set_xlim(3., 18)
    ax41.set_xlim(3., 18)
    ax41.set_xticks([4, 8, 12, 16])
    ax41.set_xticklabels(["$4$", "$8$", "$12$", "$16$"], fontsize = fontsize)

    ax12.set_xlim(3., 29)
    ax22.set_xlim(3., 29)
    ax32.set_xlim(3., 29)
    ax42.set_xlim(3., 29)
    ax42.set_xticks([5, 10, 15, 20, 25])
    ax42.set_xticklabels(["$5$", "$10$", "$15$", "$20$", "$25$"], fontsize = fontsize)

    ax12.set_ylim(0.108, 0.144)
    ax12.set_yticks([0.12, 0.14])
    ax12.set_yticklabels(["$0.12$", "$0.14$"], fontsize = fontsize)

    ax22.set_ylim(-0.15, 1.2)
    ax22.set_yticks([0., 1.])
    ax22.set_yticklabels(["$0$", "$1$"], fontsize = fontsize)
    #
    ax32.set_ylim(-0.15, 1.2)
    ax32.set_yticks([0., 1.])
    ax32.set_yticklabels(["$0$", "$1$"], fontsize = fontsize)
    #
    ax42.set_ylim(-0.175, 0.175)
    ax42.set_yticks([-0.1, 0., 0.1])
    ax42.set_yticklabels(["$-0.1$", "$0$", "$0.1$"], fontsize = fontsize)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)

    arrow = patches.FancyArrowPatch((8., .111), (8. + 1. / 0.05, .111), arrowstyle='<->', mutation_scale=8, zorder=100, linewidth=.8, color='black')
    ax12.add_patch(arrow)
    ax12.text(19., .1125, r"$ \pi / \omega_{\rm P}$", fontsize=10, color=color1, alpha=1., bbox=boxProps)

    arrow = patches.FancyArrowPatch((8., 0.118), (8. + 1. / 0.1, .118), arrowstyle='<->', mutation_scale=8, zorder=100, linewidth=.8, color='black')
    ax12.add_patch(arrow)
    ax12.text(12., .1195, r"$\pi / \omega_{\rm P}$", fontsize=10, color=color3, alpha=1., bbox=boxProps)


    legend1 = ax11.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(-.01, 1.05), edgecolor='black', ncol=1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend1 = ax12.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(-.01, 1.05), edgecolor='black', ncol=1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    # legend4 = ax4.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(0.7, 1.1), edgecolor='black', ncol=1)
    # legend4.get_frame().set_alpha(0.)
    # legend4.get_frame().set_boxstyle('Square', pad=0.1)
    # legend4.get_frame().set_linewidth(0.0)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax11.text(6., 0.147, r"$\rm Ultra$ $\rm Strong$ $\rm Coupling$", fontsize=10, color='black', alpha=1., bbox=boxProps)
    ax11.text(12., 0.1375, r"$\omega_{\rm P} = 0.5 \, \omega_{\rm phot}$", fontsize=8, color='black', alpha=1., bbox=boxProps)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax12.text(12., 0.147, r"$\rm Strong$ $\rm Coupling$", fontsize=10, color='black', alpha=1., bbox=boxProps)

    ax12.text(3., 0.147, r"$\rm{\textbf{a.)}}$", fontsize=10, color='black', alpha=1., bbox=boxProps)
    ax11.text(3., 0.147, r"$\rm{\textbf{b.)}}$", fontsize=10, color='black', alpha=1., bbox=boxProps)

    ax12.text(-2.5, 0.147, r"$X^2n_{\uparrow}n_{\downarrow}$", fontsize=12, color='black', alpha=1., bbox=boxProps)

    plt.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    #plt.show()

    plt.savefig('SCvsUSCQuad.png', format='png', bbox_inches='tight', dpi = 600)


def plotLinUSCpSC():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../data/tEvol2PhGPH50WP71WD200FD300NB10TS80WPh1230.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol2PhGPH50WP71WD156FD300NB10TS80WPh1230.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol2PhGPH50WP71WD256FD300NB10TS80WPh1230.hdf5", 'r')

    fileWP4 = h5py.File("../data/tEvol2PhGPH50WP4WD200FD300NB10TS80WPh1230.hdf5", 'r')
    fileWP5 = h5py.File("../data/tEvol2PhGPH50WP7WD200FD300NB10TS80WPh1230.hdf5", 'r')

    # readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP1 = (fileWP1['wPh'][()])[0]
    wPtWP1 = (fileWP1['wPt'][()])[0]
    wPWP1 = (fileWP1['wP'][()])[0]

    wPhWP2 = (fileWP2['wPh'][()])[0]
    wPtWP2 = (fileWP2['wPt'][()])[0]
    wPWP2 = (fileWP2['wP'][()])[0]

    wPhWP3 = (fileWP3['wPh'][()])[0]
    wPtWP3 = (fileWP3['wPt'][()])[0]
    wPWP3 = (fileWP3['wP'][()])[0]

    wPhWP4 = (fileWP4['wPh'][()])[0]
    wPtWP4 = (fileWP4['wPt'][()])[0]
    wPWP4 = (fileWP4['wP'][()])[0]

    wPhWP5 = (fileWP5['wPh'][()])[0]
    wPtWP5 = (fileWP5['wPt'][()])[0]
    wPWP5 = (fileWP5['wP'][()])[0]

    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP1, wPtWP1, wPWP1)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP1, wPtWP1, wPWP1)

    wDriveWP1 = (fileWP1['wDrive'][()])[0]
    wDriveWP2 = (fileWP2['wDrive'][()])[0]
    wDriveWP3 = (fileWP3['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDriveWP1))

    timesWP1 = fileWP1['times'][()]
    # timesWP1 = timesWP1 / (2. * np.pi / wDriveWP1)
    timesWP1 = timesWP1 / (2. * np.pi / 2.)
    pumpWP1 = fileWP1['pump'][()]

    timesWP2 = fileWP2['times'][()]
    # timesWP2 = timesWP2 / (2. * np.pi / wDriveWP2)
    timesWP2 = timesWP2 / (2. * np.pi / 2.)
    pumpWP2 = fileWP2['pump'][()]

    timesWP3 = fileWP3['times'][()]
    # timesWP3 = timesWP3 / (2. * np.pi / wDriveWP3)
    timesWP3 = timesWP3 / (2. * np.pi / 2.)
    pumpWP3 = fileWP3['pump'][()]

    timesWP4 = fileWP4['times'][()]
    timesWP4 = timesWP4 / (2. * np.pi / 2.)
    pumpWP4 = fileWP4['pump'][()]

    timesWP5 = fileWP5['times'][()]
    timesWP5 = timesWP5 / (2. * np.pi / 2.)
    pumpWP5 = fileWP5['pump'][()]

    dOccWP1 = fileWP1['dOcc'][()]
    dOccWP2 = fileWP2['dOcc'][()]
    dOccWP3 = fileWP3['dOcc'][()]
    dOccWP4 = fileWP4['dOcc'][()]
    dOccWP5 = fileWP5['dOcc'][()]

    NptWP1 = fileWP1['Npt'][()]
    NptWP2 = fileWP2['Npt'][()]
    NptWP3 = fileWP3['Npt'][()]
    NptWP4 = fileWP4['Npt'][()]
    NptWP5 = fileWP5['Npt'][()]

    NphWP1 = fileWP1['N1ph'][()]
    NphWP2 = fileWP2['N1ph'][()]
    NphWP3 = fileWP3['N1ph'][()]
    NphWP4 = fileWP4['N1ph'][()]
    NphWP5 = fileWP5['N1ph'][()]

    fig = plt.figure()
    fig.set_size_inches(6., 3.)

    gs = gridspec.GridSpec(4, 2, height_ratios=[2, 1, 1, 1])
    ax11 = fig.add_subplot(gs[0, 1])
    ax21 = fig.add_subplot(gs[1, 1], sharex=ax11)
    ax31 = fig.add_subplot(gs[2, 1], sharex=ax11)
    ax41 = fig.add_subplot(gs[3, 1], sharex=ax11)

    ax12 = fig.add_subplot(gs[0, 0], sharey=ax11)
    ax22 = fig.add_subplot(gs[1, 0], sharex=ax12, sharey=ax21)
    ax32 = fig.add_subplot(gs[2, 0], sharex=ax12, sharey=ax31)
    ax42 = fig.add_subplot(gs[3, 0], sharex=ax12, sharey=ax41)

    ax11.tick_params(direction='inout', length=4, width=.8)
    ax21.tick_params(direction='inout', length=4, width=.8)
    ax31.tick_params(direction='inout', length=4, width=.8)
    ax41.tick_params(direction='inout', length=4, width=.8)

    ax12.tick_params(direction='inout', length=4, width=.8)
    ax22.tick_params(direction='inout', length=4, width=.8)
    ax32.tick_params(direction='inout', length=4, width=.8)
    ax42.tick_params(direction='inout', length=4, width=.8)

    ax11.yaxis.tick_right()
    ax21.yaxis.tick_right()
    ax31.yaxis.tick_right()
    ax41.yaxis.tick_right()

    linewidth = 0.8

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax11.plot(timesWP1, dOccWP1 + 0.5, color=color1, label=r'$\omega_{\rm Drive} = \omega_{\rm phot}$',
              linewidth=linewidth, zorder=10)
    ax11.plot(timesWP2, dOccWP2 + 0.5, color=color2, label=r'$\omega_{\rm Drive} = \omega_-$', linewidth=linewidth,
              zorder=5)
    ax11.plot(timesWP3, dOccWP3 + 0.5, color=color3, label=r'$\omega_{\rm Drive} = \omega_+$', linewidth=linewidth)

    ax21.plot(timesWP1, 2. * NphWP1, color=color1, linewidth=linewidth, zorder=10)
    ax21.plot(timesWP2, 2. * NphWP2, color=color2, linewidth=linewidth, zorder=5)
    ax21.plot(timesWP3, 2. * NphWP3, color=color3, linewidth=linewidth)

    ax31.plot(timesWP1, NptWP1, color=color1, linewidth=linewidth, zorder=10)
    ax31.plot(timesWP2, NptWP2, color=color2, linewidth=linewidth, zorder=5)
    ax31.plot(timesWP3, NptWP3, color=color3, linewidth=linewidth)

    ax41.plot(timesWP1, pumpWP1 * 2. / 3., color=color1, label='pump', linewidth=1., zorder=10)
    ax41.plot(timesWP2, pumpWP2 * 2. / 3., color=color2, label='pump', linewidth=1., zorder=5)
    ax41.plot(timesWP3, pumpWP3 * 2. / 3., color=color3, label='pump', linewidth=1.)

    ax12.plot(timesWP4, dOccWP4 + 0.5, color=color1, label=r'$\omega_{\rm P} = \frac{\omega_{\rm phot}}{40}$',
              linewidth=linewidth)
    ax12.plot(timesWP5, dOccWP5 + 0.5, color=color3, label=r'$\omega_{\rm P} = \frac{\omega_{\rm phot}}{20}$',
              linewidth=linewidth)

    ax22.plot(timesWP4, 2. * NphWP4, color=color1, linewidth=linewidth)
    ax22.plot(timesWP5, 2. * NphWP5, color=color3, linewidth=linewidth)

    ax32.plot(timesWP4, NptWP4, color=color1, linewidth=linewidth)
    ax32.plot(timesWP5, NptWP5, color=color3, linewidth=linewidth)

    ax42.plot(timesWP4, pumpWP4 * 2. / 3., color=color1, label='pump', linewidth=1.)
    ax42.plot(timesWP5, pumpWP5 * 2. / 3., color=color3, label='pump', linewidth=1., linestyle='--', dashes=[4, 4])

    ax12.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize, labelpad = 3.)
    ax22.set_ylabel(r"$N_{\rm phon}$", fontsize=fontsize, labelpad = 15.)
    ax32.set_ylabel(r"$N_{\rm phot}$", fontsize=fontsize, labelpad = 15.)
    ax42.set_ylabel(r"$F(t) / F_0$", fontsize=fontsize, labelpad = 0.)

    ax41.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm phot}]$", fontsize=fontsize)
    ax42.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm phot}]$", fontsize=fontsize)

    plt.setp(ax11.get_xticklabels(), visible=False)
    plt.setp(ax21.get_xticklabels(), visible=False)
    plt.setp(ax31.get_xticklabels(), visible=False)

    plt.setp(ax12.get_xticklabels(), visible=False)
    plt.setp(ax22.get_xticklabels(), visible=False)
    plt.setp(ax32.get_xticklabels(), visible=False)

    plt.setp(ax11.get_yticklabels(), visible=False)
    plt.setp(ax21.get_yticklabels(), visible=False)
    plt.setp(ax31.get_yticklabels(), visible=False)
    plt.setp(ax41.get_yticklabels(), visible=False)

    ax11.set_xlim(3., 18)
    ax21.set_xlim(3., 18)
    ax31.set_xlim(3., 18)
    ax41.set_xlim(3., 18)
    ax41.set_xticks([4, 8, 12, 16])
    ax41.set_xticklabels(["$4$", "$8$", "$12$", "$16$"], fontsize = fontsize)

    ax12.set_xlim(3., 29)
    ax22.set_xlim(3., 29)
    ax32.set_xlim(3., 29)
    ax42.set_xlim(3., 29)
    ax42.set_xticks([5, 10, 15, 20, 25])
    ax42.set_xticklabels(["$5$", "$10$", "$15$", "$20$", "$25$"], fontsize = fontsize)

    # ax1.set_ylim(0.08, 0.21)
    ax12.set_yticks([0.11, 0.12])
    ax12.set_yticklabels(["$0.11$", "$0.12$"], fontsize = fontsize)
    #
    ax22.set_ylim(-0.15, 1.7)
    ax22.set_yticks([0, 1])
    ax22.set_yticklabels(["$0$", "$1$", '$2$'], fontsize = fontsize)
    #
    ax32.set_ylim(-0.15, 1.4)
    ax32.set_yticks([0, 1])
    ax32.set_yticklabels(["$0$", "$1$"], fontsize = fontsize)
    #
    ax42.set_ylim(-0.175, 0.175)
    ax42.set_yticks([-0.1, 0., 0.1])
    ax42.set_yticklabels(["$-0.1$", "$0$", "$0.1$"], fontsize = fontsize)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)

    arrow = patches.FancyArrowPatch((7., .11), (7. + 1. / 0.05, .11), arrowstyle='<->', mutation_scale=8, zorder=100, linewidth=.8, color='black')
    ax12.add_patch(arrow)
    ax12.text(18., .1106, r"$\pi/\omega_{\rm P}$", fontsize=10, color=color1, alpha=1., bbox=boxProps)

    arrow = patches.FancyArrowPatch((7., 0.113), (7. + 1. / 0.1, .113), arrowstyle='<->', mutation_scale=8, zorder=100, linewidth=.8, color='black')
    ax12.add_patch(arrow)
    ax12.text(12., .1138, r"$\pi / \omega_{\rm P}$", fontsize=10, color=color3, alpha=1., bbox=boxProps)

    legend1 = ax11.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(-.01, 1.05), edgecolor='black',
                          ncol=1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend1 = ax12.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(-.01, 1.05), edgecolor='black',
                          ncol=1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    # legend4 = ax4.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(0.7, 1.1), edgecolor='black', ncol=1)
    # legend4.get_frame().set_alpha(0.)
    # legend4.get_frame().set_boxstyle('Square', pad=0.1)
    # legend4.get_frame().set_linewidth(0.0)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax11.text(7., 0.124, r"$\rm Ultra$ $\rm Strong$ $\rm Coupling$", fontsize=10, color='black', alpha=1., bbox=boxProps)
    ax11.text(12., 0.12, r"$\omega_{\rm P} = 0.5 \, \omega_{\rm phot}$", fontsize=8, color='black', alpha=1.,
              bbox=boxProps)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax12.text(12., 0.124, r"$\rm Strong$ $\rm Coupling$", fontsize=10, color='black', alpha=1., bbox=boxProps)


    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax12.text(3., 0.124, r"$\rm{\textbf{a.)}}$", fontsize=10, color='black', alpha=1., bbox=boxProps)
    ax11.text(3., 0.124, r"$\rm{\textbf{b.)}}$", fontsize=10, color='black', alpha=1., bbox=boxProps)

    ax12.text(-2.5, 0.124, r"$X^2n$", fontsize=12, color='black', alpha=1., bbox=boxProps)

    plt.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    #plt.show()

    plt.savefig('SCvsUSCLin.png', format='png', bbox_inches='tight', dpi=600)

def plotAsOfWPLin():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../data/tEvol2PhGPH50WP0WD200FD200NB10TS80.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol2PhGPH50WP7WD200FD200NB10TS80.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD200NB10TS80.hdf5", 'r')
    fileWP4 = h5py.File("../data/tEvol2PhGPH50WP28WD200FD200NB10TS80.hdf5", 'r')

    ## read in stuff
    #fileWP1 = h5py.File("../data/tEvol2PhGPH50WP0WD200FD200NB10TS80.hdf5", 'r')
    #fileWP2 = h5py.File("../data/tEvol2PhGPH50WP10WD200FD200NB10TS80.hdf5", 'r')
    #fileWP3 = h5py.File("../data/tEvol2PhGPH50WP20WD200FD200NB10TS80.hdf5", 'r')
    #fileWP4 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB10TS80.hdf5", 'r')

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
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

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
    ax1.set_xticks([4, 8, 12, 16])
    ax2.set_xticks([4, 8, 12, 16])
    ax3.set_xticks([4, 8, 12, 16])
    ax3.set_xticklabels(["$4$", "$8$", "$12$", "$16$"], fontsize = fontsize)

    # ax1.set_ylim(0.08, 0.21)
    ax1.set_yticks([0.11, 0.12, 0.13])
    ax1.set_yticklabels(["$0.11$", "$0.12$", "$0.13$"], fontsize = fontsize)
    #
    ax2.set_yticks([0., 0.3, 0.6])
    ax2.set_yticklabels(["$0$", "$0.3$", '$0.6$'], fontsize = fontsize)
    #
    ax3.set_yticks([0, 1, 2])
    ax3.set_yticklabels(["$0$", "$1$", "$2$"], fontsize = fontsize)
    #
    ax4.set_ylim(-0.4, 0.4)
    ax4.set_yticks([-0.2, 0., 0.2])
    ax4.set_yticklabels(["$-0.2$", "$0$", "$0.2$"], fontsize = fontsize)


    legend1 = ax1.legend(fontsize = fontsize - 4, loc = 'upper left', bbox_to_anchor=(.0, 1.), edgecolor = 'black', ncol = 1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(0.7, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    #plt.savefig('functionOfWPLin.png', format='png', bbox_inches='tight', dpi = 600)

def plot1Phonvs2Phon():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileC = h5py.File("../data/tEvol1PhGPH10WP20WD200FD300NB10TS80WPh2011.hdf5", 'r')
    fileQ = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB10TS80WPh2022.hdf5", 'r')

    # read in stuff
    #fileC = h5py.File("../data/tEvol1PhGPH10WP71WD200FD200NB10.hdf5", 'r')
    #fileQ = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB10TS80.hdf5", 'r')

    # readInPrmsAndAssert(fileC, fileQ1)
    wPhQ = (fileQ['wPh'][()])[0]
    wPtQ = (fileQ['wPt'][()])[0]
    wPQ = (fileQ['wP'][()])[0]

    wDriveC = (fileC['wDrive'][()])[0]
    wDriveQ = (fileQ['wDrive'][()])[0]


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
    dOccQ = fileQ['dOcc'][()]

    Npt = fileC['Npt'][()]
    Nph = fileC['N1ph'][()]

    NptQ = fileQ['Npt'][()]
    NphQ = fileQ['N1ph'][()]

    fig = plt.figure()
    fig.set_size_inches(5., 3.)

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex = ax1)
    ax3 = fig.add_subplot(gs[2], sharex = ax1)
    ax4 = fig.add_subplot(gs[3], sharex = ax1)

    ax1.tick_params(direction='inout', length=4, width=.8)
    ax2.tick_params(direction='inout', length=4, width=.8)
    ax3.tick_params(direction='inout', length=4, width=.8)
    ax4.tick_params(direction='inout', length=4, width=.8)

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'

    ax1.plot(times, dOccQ + 0.5, color=color3, label='2Phonons', linewidth = 1.)
    ax1.plot(times, dOccC + 0.5, color=color1, label='1Phonon', linewidth = 1., zorder = 10)

    ax2.plot(times, NphQ, color=color3, linewidth = 1.)
    ax2.plot(times, Nph, color=color1, linewidth = 1., zorder = 10)

    ax3.plot(times, NptQ, color=color3, label='Driven Cavity', linewidth = 1.)
    ax3.plot(times, Npt, color=color1, label='Classical Driving', linewidth = 1., zorder = 10)

    ax4.plot(times, pump * 2 / 3, color='red', label='pump', linewidth = 1.)


    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize, labelpad=3.)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize = fontsize, labelpad = 7.)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize = fontsize, labelpad = 7.)
    ax4.set_ylabel(r"$F(t)/F_0$", fontsize = fontsize, labelpad=0.)


    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm phot}]$", fontsize = fontsize)


    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)

    ax1.set_xlim(3., 23)
    ax2.set_xlim(3., 23)
    ax3.set_xlim(3., 23)
    ax4.set_xlim(3., 23)
    ax1.set_xticks([5, 10, 15, 20])
    ax2.set_xticks([5, 10, 15, 20])
    ax3.set_xticks([5, 10, 15, 20])
    ax4.set_xticks([5, 10, 15, 20])
    ax4.set_xticklabels(["$5$", "$10$", "$15$", "$20$"], fontsize = fontsize)


    ax1.set_ylim(0.105, 0.135)
    ax1.set_yticks([0.11, 0.12, 0.13])
    ax1.set_yticklabels(["$0.11$", "$0.12$", "$0.13$"], fontsize = fontsize)

    ax2.set_ylim(-0.1, 0.7)
    ax2.set_yticks([0, 0.5])
    ax2.set_yticklabels(["$0$", "$0.5$"], fontsize = fontsize)
    ax3.set_ylim(-0.1, 0.9)
    ax3.set_yticks([0, 0.5])
    ax3.set_yticklabels(["$0$", "$0.5$"], fontsize = fontsize)
    ax4.set_ylim(-0.175, 0.175)
    ax4.set_yticks([-0.1, 0., 0.1])
    ax4.set_yticklabels(["$-0.1$", "$0$", "$0.1$"], fontsize = fontsize)


    legend1 = ax1.legend(fontsize = fontsize - 2, loc = 'upper left', bbox_to_anchor=(.0, 1.), edgecolor = 'black', ncol = 1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    # legend2 = ax2.legend(fontsize = fontsize - 4, loc = 'upper left', bbox_to_anchor=(.45, 1.), edgecolor = 'black', ncol = 1)
    # legend2.get_frame().set_alpha(0.)
    # legend2.get_frame().set_boxstyle('Square', pad=0.1)
    # legend2.get_frame().set_linewidth(0.0)


    legend4 = ax4.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(.45, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    #plt.show()

    plt.savefig('1PhononVS2Phonon.png', format='png', bbox_inches='tight', dpi = 600)




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