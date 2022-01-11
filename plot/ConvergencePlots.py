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

def plotAsOfWPLinDTConv():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../data/tEvol2PhGPH50WP20WD200FD200NB6TS20.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol2PhGPH25WP20WD200FD200NB6TS20.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol2PhGPH50WP20WD200FD200NB6TS20.hdf5", 'r')
    fileWP4 = h5py.File("../data/tEvol2PhGPH50WP20WD200FD200NB6TS20.hdf5", 'r')
    fileWP5 = h5py.File("../data/tEvol2PhGPH50WP20WD200FD200NB6TS20.hdf5", 'r')

    # readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP4 = (fileWP3['wPh'][()])[0]
    wPtWP4 = (fileWP3['wPt'][()])[0]
    wPWP4 = (fileWP3['wP'][()])[0]




    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP4, wPtWP4, wPWP4)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP4, wPtWP4, wPWP4)

    wDrive = (fileWP1['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDrive))


    times1 = fileWP1['times'][()]
    times1 = times1 / (2. * np.pi / wDrive)
    pump1 = fileWP1['pump'][()]

    times2 = fileWP2['times'][()]
    times2 = times2 / (2. * np.pi / wDrive)
    pump2 = fileWP2['pump'][()]

    times3 = fileWP3['times'][()]
    times3 = times3 / (2. * np.pi / wDrive)
    pump3 = fileWP3['pump'][()]

    times4 = fileWP4['times'][()]
    times4 = times4 / (2. * np.pi / wDrive)
    pump4 = fileWP4['pump'][()]

    times5 = fileWP5['times'][()]
    times5 = times5 / (2. * np.pi / wDrive)
    pump5 = fileWP5['pump'][()]

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
    fig.set_size_inches(3., 3.)
    # fig.set_size_inches(5., 5.)

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex = ax1)
    ax3 = fig.add_subplot(gs[2], sharex = ax1)
    ax4 = fig.add_subplot(gs[3], sharex = ax1)

    ax1.tick_params(direction='inout', length=6, width=.5)
    ax2.tick_params(direction='inout', length=6, width=.5)
    ax3.tick_params(direction='inout', length=6, width=.5)

    linewidth = 1.2

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax1.plot(times1, dOccWP1 + 0.5, color=color1, label=r'$N_t = 10$', linewidth = linewidth, zorder = 2)
    ax1.plot(times2, dOccWP2 + 0.5, color=color2, label=r'$N_t = 20$', linewidth = linewidth, zorder = 10)
    ax1.plot(times3, dOccWP3 + 0.5, color=color3, label=r'$N_t = 40$', linewidth = linewidth, zorder = 4)
    ax1.plot(times4, dOccWP4 + 0.5, color=color4, label=r'$N_t = 120$', linewidth = linewidth, zorder = 5)
    ax1.plot(times5, dOccWP5 + 0.5, color=color5, label=r'$N_t = 80$', linewidth = linewidth, zorder = 5, linestyle = '--', dashes=[4, 4])

    ax2.plot(times1, NphWP1, color=color1, linewidth = linewidth)
    ax2.plot(times2, NphWP2, color=color2, linewidth = linewidth)
    ax2.plot(times3, NphWP3, color=color3, linewidth = linewidth)
    ax2.plot(times4, NphWP4, color=color4, linewidth = linewidth)
    ax2.plot(times5, NphWP5, color=color5, linewidth = linewidth, linestyle = '--', dashes=[4, 4])

    ax3.plot(times1, NptWP1, color=color1, linewidth = linewidth)
    ax3.plot(times2, NptWP2, color=color2, linewidth = linewidth)
    ax3.plot(times3, NptWP3, color=color3, linewidth = linewidth)
    ax3.plot(times4, NptWP4, color=color4, linewidth = linewidth)
    ax3.plot(times5, NptWP5, color=color5, linewidth = linewidth, linestyle = '--', dashes=[4, 4])

    # ax4.plot(times1, pump1, color=colorPump, label='pump', linewidth = 1.)
    # ax4.plot(times2, pump2, color=colorPump, label='pump', linewidth = 1.)
    ax4.plot(times4, pump4, color=colorPump, label='pump', linewidth = 1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize = fontsize)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize = fontsize)
    ax4.set_ylabel(r"$F(t)$", fontsize = fontsize)


    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize = fontsize)


    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    # ax1.set_xlim(10., 14.)
    # ax2.set_xlim(10., 14.)
    # ax3.set_xlim(10., 14.)
    # ax4.set_xlim(10., 14.)

    # ax1.set_ylim(0.11, 0.115)
    # ax2.set_ylim(0.05, 0.25)
    # ax3.set_ylim(0.08, 0.22)
    # ax4.set_ylim(-0.25, 0.25)
    # ax1.set_xticks([4, 6, 8, 10, 12])
    # ax2.set_xticks([4, 6, 8, 10, 12])
    # ax3.set_xticks([4, 6, 8, 10, 12])
    # ax3.set_xticklabels(["$4$", "$6$", "$8$", "$10$", "$12$"])

    # ax1.set_ylim(0.08, 0.21)
    # ax1.set_yticks([0.1, 0.13, 0.16])
    # ax1.set_yticklabels(["$0.1$", "$0.13$", "$0.16$"])
    #
    # ax2.set_yticks([0, 1])
    # ax2.set_yticklabels(["$0$", "$1$"])
    #
    # ax3.set_yticks([0, 1, 2])
    # ax3.set_yticklabels(["$0$", "$1$", "$2$"])
    #
    ax4.set_yticks([-0.2, 0., 0.2])
    ax4.set_yticklabels(["$-0.2$", "$0$", "$0.2$"])


    legend1 = ax1.legend(fontsize = fontsize - 4, loc = 'upper left', bbox_to_anchor=(.0, 1.), edgecolor = 'black', ncol = 3)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize, loc='upper left', bbox_to_anchor=(0.6, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    # plt.savefig('functionOfWPLinConv.png', format='png', bbox_inches='tight', dpi = 600)


def plotAsOfWPQuadDTConv():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB6TS10.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB6TS20.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB6TS40.hdf5", 'r')
    fileWP4 = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB6TS120.hdf5", 'r')
    fileWP5 = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB6TS80.hdf5", 'r')

    # readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP4 = (fileWP3['wPh'][()])[0]
    wPtWP4 = (fileWP3['wPt'][()])[0]
    wPWP4 = (fileWP3['wP'][()])[0]

    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP4, wPtWP4, wPWP4)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP4, wPtWP4, wPWP4)

    wDrive = (fileWP1['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDrive))

    times1 = fileWP1['times'][()]
    times1 = times1 / (2. * np.pi / wDrive)
    pump1 = fileWP1['pump'][()]

    times2 = fileWP2['times'][()]
    times2 = times2 / (2. * np.pi / wDrive)
    pump2 = fileWP2['pump'][()]

    times3 = fileWP3['times'][()]
    times3 = times3 / (2. * np.pi / wDrive)
    pump3 = fileWP3['pump'][()]

    times4 = fileWP4['times'][()]
    times4 = times4 / (2. * np.pi / wDrive)
    pump4 = fileWP4['pump'][()]

    times5 = fileWP5['times'][()]
    times5 = times5 / (2. * np.pi / wDrive)
    pump5 = fileWP5['pump'][()]

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
    fig.set_size_inches(3., 3.)
    # fig.set_size_inches(5., 5.)

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    ax4 = fig.add_subplot(gs[3], sharex=ax1)

    ax1.tick_params(direction='inout', length=6, width=.5)
    ax2.tick_params(direction='inout', length=6, width=.5)
    ax3.tick_params(direction='inout', length=6, width=.5)

    linewidth = 1.2

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax1.plot(times1, dOccWP1 + 0.5, color=color1, label=r'$N_t = 10$', linewidth=linewidth, zorder=2)
    ax1.plot(times2, dOccWP2 + 0.5, color=color2, label=r'$N_t = 20$', linewidth=linewidth, zorder=3)
    ax1.plot(times3, dOccWP3 + 0.5, color=color3, label=r'$N_t = 40$', linewidth=linewidth, zorder=4)
    ax1.plot(times4, dOccWP4 + 0.5, color=color4, label=r'$N_t = 120$', linewidth=linewidth, zorder=5)
    ax1.plot(times5, dOccWP5 + 0.5, color=color5, label=r'$N_t = 80$', linewidth=linewidth, zorder=5, linestyle='--',
             dashes=[4, 4])

    ax2.plot(times1, NphWP1, color=color1, linewidth=linewidth)
    ax2.plot(times2, NphWP2, color=color2, linewidth=linewidth)
    ax2.plot(times3, NphWP3, color=color3, linewidth=linewidth)
    ax2.plot(times4, NphWP4, color=color4, linewidth=linewidth)
    ax2.plot(times5, NphWP5, color=color5, linewidth=linewidth, linestyle='--', dashes=[4, 4])

    ax3.plot(times1, NptWP1, color=color1, linewidth=linewidth)
    ax3.plot(times2, NptWP2, color=color2, linewidth=linewidth)
    ax3.plot(times3, NptWP3, color=color3, linewidth=linewidth)
    ax3.plot(times4, NptWP4, color=color4, linewidth=linewidth)
    ax3.plot(times5, NptWP5, color=color5, linewidth=linewidth, linestyle='--', dashes=[4, 4])

    # ax4.plot(times1, pump1, color=colorPump, label='pump', linewidth = 1.)
    # ax4.plot(times2, pump2, color=colorPump, label='pump', linewidth = 1.)
    ax4.plot(times4, pump4, color=colorPump, label='pump', linewidth=1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize=fontsize)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize=fontsize)
    ax4.set_ylabel(r"$F(t)$", fontsize=fontsize)

    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize=fontsize)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim(10., 14.)
    ax2.set_xlim(10., 14.)
    ax3.set_xlim(10., 14.)
    ax4.set_xlim(10., 14.)

    ax1.set_ylim(0.116, 0.118)
    ax2.set_ylim(0., 0.05)
    ax3.set_ylim(0., 0.07)
    ax4.set_ylim(-0.25, 0.25)

    # ax1.set_xticks([4, 6, 8, 10, 12])
    # ax2.set_xticks([4, 6, 8, 10, 12])
    # ax3.set_xticks([4, 6, 8, 10, 12])
    # ax3.set_xticklabels(["$4$", "$6$", "$8$", "$10$", "$12$"])

    # ax1.set_ylim(0.08, 0.21)
    # ax1.set_yticks([0.1, 0.13, 0.16])
    # ax1.set_yticklabels(["$0.1$", "$0.13$", "$0.16$"])
    #
    # ax2.set_yticks([0, 1])
    # ax2.set_yticklabels(["$0$", "$1$"])
    #
    # ax3.set_yticks([0, 1, 2])
    # ax3.set_yticklabels(["$0$", "$1$", "$2$"])
    #
    ax4.set_yticks([-0.2, 0., 0.2])
    ax4.set_yticklabels(["$-0.2$", "$0$", "$0.2$"])

    legend1 = ax1.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(.0, 1.), edgecolor='black', ncol=3)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(0.6, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    # plt.show()

    plt.savefig('functionOfWPQuadConv.png', format='png', bbox_inches='tight', dpi=600)


def plotAsOfWPLinDTConvNB():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB4.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB6TS80.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB8TS80.hdf5", 'r')
    fileWP4 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB10TS80.hdf5", 'r')
    fileWP5 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB10TS80.hdf5", 'r')

    # readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP4 = (fileWP3['wPh'][()])[0]
    wPtWP4 = (fileWP3['wPt'][()])[0]
    wPWP4 = (fileWP3['wP'][()])[0]

    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP4, wPtWP4, wPWP4)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP4, wPtWP4, wPWP4)

    wDrive = (fileWP1['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDrive))

    times1 = fileWP1['times'][()]
    times1 = times1 / (2. * np.pi / wDrive)
    pump1 = fileWP1['pump'][()]

    times2 = fileWP2['times'][()]
    times2 = times2 / (2. * np.pi / wDrive)
    pump2 = fileWP2['pump'][()]

    times3 = fileWP3['times'][()]
    times3 = times3 / (2. * np.pi / wDrive)
    pump3 = fileWP3['pump'][()]

    times4 = fileWP4['times'][()]
    times4 = times4 / (2. * np.pi / wDrive)
    pump4 = fileWP4['pump'][()]

    times5 = fileWP5['times'][()]
    times5 = times5 / (2. * np.pi / wDrive)
    pump5 = fileWP5['pump'][()]

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
    fig.set_size_inches(3., 3.)
    # fig.set_size_inches(5., 5.)

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    ax4 = fig.add_subplot(gs[3], sharex=ax1)

    ax1.tick_params(direction='inout', length=6, width=.5)
    ax2.tick_params(direction='inout', length=6, width=.5)
    ax3.tick_params(direction='inout', length=6, width=.5)

    linewidth = 1.2

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax1.plot(times1, dOccWP1 + 0.5, color=color1, label=r'$N_t = 10$', linewidth=linewidth, zorder=2)
    ax1.plot(times2, dOccWP2 + 0.5, color=color2, label=r'$N_t = 20$', linewidth=linewidth, zorder=3)
    ax1.plot(times3, dOccWP3 + 0.5, color=color3, label=r'$N_t = 40$', linewidth=linewidth, zorder=4)
    ax1.plot(times4, dOccWP4 + 0.5, color=color4, label=r'$N_t = 120$', linewidth=linewidth, zorder=5)
    ax1.plot(times5, dOccWP5 + 0.5, color=color5, label=r'$N_t = 80$', linewidth=linewidth, zorder=5, linestyle='--',
             dashes=[4, 4])

    ax2.plot(times1, NphWP1, color=color1, linewidth=linewidth)
    ax2.plot(times2, NphWP2, color=color2, linewidth=linewidth)
    ax2.plot(times3, NphWP3, color=color3, linewidth=linewidth)
    ax2.plot(times4, NphWP4, color=color4, linewidth=linewidth)
    ax2.plot(times5, NphWP5, color=color5, linewidth=linewidth, linestyle='--', dashes=[4, 4])

    ax3.plot(times1, NptWP1, color=color1, linewidth=linewidth)
    ax3.plot(times2, NptWP2, color=color2, linewidth=linewidth)
    ax3.plot(times3, NptWP3, color=color3, linewidth=linewidth)
    ax3.plot(times4, NptWP4, color=color4, linewidth=linewidth)
    ax3.plot(times5, NptWP5, color=color5, linewidth=linewidth, linestyle='--', dashes=[4, 4])

    # ax4.plot(times1, pump1, color=colorPump, label='pump', linewidth = 1.)
    # ax4.plot(times2, pump2, color=colorPump, label='pump', linewidth = 1.)
    ax4.plot(times4, pump4, color=colorPump, label='pump', linewidth=1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize=fontsize)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize=fontsize)
    ax4.set_ylabel(r"$F(t)$", fontsize=fontsize)

    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize=fontsize)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim(0., 20.)
    ax2.set_xlim(0., 20.)
    ax3.set_xlim(0., 20.)
    ax4.set_xlim(0., 20.)

    # ax1.set_ylim(0.11, 0.115)
    # ax2.set_ylim(0.05, 0.25)
    # ax3.set_ylim(0.08, 0.22)
    # ax4.set_ylim(-0.25, 0.25)
    # ax1.set_xticks([4, 6, 8, 10, 12])
    # ax2.set_xticks([4, 6, 8, 10, 12])
    # ax3.set_xticks([4, 6, 8, 10, 12])
    # ax3.set_xticklabels(["$4$", "$6$", "$8$", "$10$", "$12$"])

    # ax1.set_ylim(0.08, 0.21)
    # ax1.set_yticks([0.1, 0.13, 0.16])
    # ax1.set_yticklabels(["$0.1$", "$0.13$", "$0.16$"])
    #
    # ax2.set_yticks([0, 1])
    # ax2.set_yticklabels(["$0$", "$1$"])
    #
    # ax3.set_yticks([0, 1, 2])
    # ax3.set_yticklabels(["$0$", "$1$", "$2$"])
    #
    ax4.set_yticks([-0.2, 0., 0.2])
    ax4.set_yticklabels(["$-0.2$", "$0$", "$0.2$"])

    legend1 = ax1.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(.0, 1.), edgecolor='black', ncol=3)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(0.6, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    # plt.savefig('functionOfWPLinConv.png', format='png', bbox_inches='tight', dpi=600)


def plotAsOfWPQuadDTConvNB():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB4.hdf5", 'r')
    fileWP2 = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB6TS80.hdf5", 'r')
    fileWP3 = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB8TS80.hdf5", 'r')
    fileWP4 = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB10TS80.hdf5", 'r')
    fileWP5 = h5py.File("../data/tEvol2PhGPH20WP50WD200FD200NB10TS80.hdf5", 'r')

    # readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP4 = (fileWP3['wPh'][()])[0]
    wPtWP4 = (fileWP3['wPt'][()])[0]
    wPWP4 = (fileWP3['wP'][()])[0]

    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP4, wPtWP4, wPWP4)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP4, wPtWP4, wPWP4)

    wDrive = (fileWP1['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDrive))

    times1 = fileWP1['times'][()]
    times1 = times1 / (2. * np.pi / wDrive)
    pump1 = fileWP1['pump'][()]

    times2 = fileWP2['times'][()]
    times2 = times2 / (2. * np.pi / wDrive)
    pump2 = fileWP2['pump'][()]

    times3 = fileWP3['times'][()]
    times3 = times3 / (2. * np.pi / wDrive)
    pump3 = fileWP3['pump'][()]

    times4 = fileWP4['times'][()]
    times4 = times4 / (2. * np.pi / wDrive)
    pump4 = fileWP4['pump'][()]

    times5 = fileWP5['times'][()]
    times5 = times5 / (2. * np.pi / wDrive)
    pump5 = fileWP5['pump'][()]

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
    fig.set_size_inches(3., 3.)
    # fig.set_size_inches(5., 5.)

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    ax4 = fig.add_subplot(gs[3], sharex=ax1)

    ax1.tick_params(direction='inout', length=6, width=.5)
    ax2.tick_params(direction='inout', length=6, width=.5)
    ax3.tick_params(direction='inout', length=6, width=.5)

    linewidth = 1.2

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax1.plot(times1, dOccWP1 + 0.5, color=color1, label=r'$N_t = 10$', linewidth=linewidth, zorder=2)
    ax1.plot(times2, dOccWP2 + 0.5, color=color2, label=r'$N_t = 20$', linewidth=linewidth, zorder=3)
    ax1.plot(times3, dOccWP3 + 0.5, color=color3, label=r'$N_t = 40$', linewidth=linewidth, zorder=4)
    ax1.plot(times4, dOccWP4 + 0.5, color=color4, label=r'$N_t = 120$', linewidth=linewidth, zorder=5)
    ax1.plot(times5, dOccWP5 + 0.5, color=color5, label=r'$N_t = 80$', linewidth=linewidth, zorder=5, linestyle='--',
             dashes=[4, 4])

    ax2.plot(times1, NphWP1, color=color1, linewidth=linewidth)
    ax2.plot(times2, NphWP2, color=color2, linewidth=linewidth)
    ax2.plot(times3, NphWP3, color=color3, linewidth=linewidth)
    ax2.plot(times4, NphWP4, color=color4, linewidth=linewidth)
    ax2.plot(times5, NphWP5, color=color5, linewidth=linewidth, linestyle='--', dashes=[4, 4])

    ax3.plot(times1, NptWP1, color=color1, linewidth=linewidth)
    ax3.plot(times2, NptWP2, color=color2, linewidth=linewidth)
    ax3.plot(times3, NptWP3, color=color3, linewidth=linewidth)
    ax3.plot(times4, NptWP4, color=color4, linewidth=linewidth)
    ax3.plot(times5, NptWP5, color=color5, linewidth=linewidth, linestyle='--', dashes=[4, 4])

    # ax4.plot(times1, pump1, color=colorPump, label='pump', linewidth = 1.)
    # ax4.plot(times2, pump2, color=colorPump, label='pump', linewidth = 1.)
    ax4.plot(times4, pump4, color=colorPump, label='pump', linewidth=1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize=fontsize)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize=fontsize)
    ax4.set_ylabel(r"$F(t)$", fontsize=fontsize)

    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize=fontsize)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim(0., 20.)
    ax2.set_xlim(0., 20.)
    ax3.set_xlim(0., 20.)
    ax4.set_xlim(0., 20.)

    # ax1.set_ylim(0.11, 0.115)
    # ax2.set_ylim(0.05, 0.25)
    # ax3.set_ylim(0.08, 0.22)
    # ax4.set_ylim(-0.25, 0.25)
    # ax1.set_xticks([4, 6, 8, 10, 12])
    # ax2.set_xticks([4, 6, 8, 10, 12])
    # ax3.set_xticks([4, 6, 8, 10, 12])
    # ax3.set_xticklabels(["$4$", "$6$", "$8$", "$10$", "$12$"])

    # ax1.set_ylim(0.08, 0.21)
    # ax1.set_yticks([0.1, 0.13, 0.16])
    # ax1.set_yticklabels(["$0.1$", "$0.13$", "$0.16$"])
    #
    # ax2.set_yticks([0, 1])
    # ax2.set_yticklabels(["$0$", "$1$"])
    #
    # ax3.set_yticks([0, 1, 2])
    # ax3.set_yticklabels(["$0$", "$1$", "$2$"])
    #
    ax4.set_yticks([-0.2, 0., 0.2])
    ax4.set_yticklabels(["$-0.2$", "$0$", "$0.2$"])

    legend1 = ax1.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(.0, 1.), edgecolor='black', ncol=3)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(0.6, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    # plt.savefig('functionOfWPLinConv.png', format='png', bbox_inches='tight', dpi=600)


def plotConvInDOccQuad():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileNB1 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB4TS40WPh2022.hdf5", 'r')
    fileNB2 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB6TS40WPh2022.hdf5", 'r')
    fileNB3 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB8TS40WPh2022.hdf5", 'r')
    fileNB4 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB10TS40WPh2022.hdf5", 'r')
    fileNB5 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB12TS40WPh2022.hdf5", 'r')

    # read in stuff
    fileNT1 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB8TS10WPh2022.hdf5", 'r')
    fileNT2 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB8TS20WPh2022.hdf5", 'r')
    fileNT3 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB8TS40WPh2022.hdf5", 'r')
    fileNT4 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB8TS80WPh2022.hdf5", 'r')
    fileNT5 = h5py.File("../data/tEvol2PhGPH20WP14WD200FD300NB8TS160WPh2022.hdf5", 'r')

    # readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP4 = (fileNB3['wPh'][()])[0]
    wPtWP4 = (fileNB3['wPt'][()])[0]
    wPWP4 = (fileNB3['wP'][()])[0]

    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP4, wPtWP4, wPWP4)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP4, wPtWP4, wPWP4)

    wDrive = (fileNB1['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDrive))

    timesNB = fileNB1['times'][()]
    timesNB = timesNB / (2. * np.pi / wDrive)
    pump1 = fileNB1['pump'][()]

    dOccNB1 = fileNB1['dOcc'][()]
    dOccNB2 = fileNB2['dOcc'][()]
    dOccNB3 = fileNB3['dOcc'][()]
    dOccNB4 = fileNB4['dOcc'][()]
    dOccNB5 = fileNB5['dOcc'][()]

    times1NT = fileNT1['times'][()]
    times1NT = times1NT / (2. * np.pi / wDrive)
    times2NT = fileNT2['times'][()]
    times2NT = times2NT / (2. * np.pi / wDrive)
    times3NT = fileNT3['times'][()]
    times3NT = times3NT / (2. * np.pi / wDrive)
    times4NT = fileNT4['times'][()]
    times4NT = times4NT / (2. * np.pi / wDrive)
    times5NT = fileNT5['times'][()]
    times5NT = times5NT / (2. * np.pi / wDrive)

    dOccNT1 = fileNT1['dOcc'][()]
    dOccNT2 = fileNT2['dOcc'][()]
    dOccNT3 = fileNT3['dOcc'][()]
    dOccNT4 = fileNT4['dOcc'][()]
    dOccNT5 = fileNT5['dOcc'][()]

    fig = plt.figure()
    fig.set_size_inches(3., 3.)
    # fig.set_size_inches(5., 5.)

    gs = gridspec.GridSpec(3, 1, height_ratios=[2, 2, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax4 = fig.add_subplot(gs[2], sharex=ax1)

    ax1.tick_params(direction='inout', length=4, width=.8)
    ax2.tick_params(direction='inout', length=4, width=.8)
    ax4.tick_params(direction='inout', length=4, width=.8)

    linewidth = 1.2

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax1.plot(timesNB, dOccNB1 + 0.5, color=color1, label=r'$N_B = 4$', linewidth=linewidth, zorder=2)
    ax1.plot(timesNB, dOccNB2 + 0.5, color=color2, label=r'$N_B = 6$', linewidth=linewidth, zorder=3)
    ax1.plot(timesNB, dOccNB3 + 0.5, color=color3, label=r'$N_B = 8$', linewidth=linewidth, zorder=4)
    ax1.plot(timesNB, dOccNB4 + 0.5, color=color4, label=r'$N_B = 10$', linewidth=linewidth, zorder=5)
    ax1.plot(timesNB, dOccNB5 + 0.5, color=color5, label=r'$N_B = 12$', linewidth=linewidth, zorder=5, linestyle='--',
             dashes=[4, 4])

    ax2.plot(times1NT, dOccNT1 + 0.5, color=color1, label=r'$N_t = 10$', linewidth=linewidth)
    ax2.plot(times2NT, dOccNT2 + 0.5, color=color2, label=r'$N_t = 20$', linewidth=linewidth)
    ax2.plot(times3NT, dOccNT3 + 0.5, color=color3, label=r'$N_t = 40$', linewidth=linewidth)
    ax2.plot(times4NT, dOccNT4 + 0.5, color=color4, label=r'$N_t = 80$', linewidth=linewidth)
    ax2.plot(times5NT, dOccNT5 + 0.5, color=color5, label=r'$N_t = 160$', linewidth=linewidth, linestyle='--', dashes=[4, 4])

    # ax4.plot(times1, pump1, color=colorPump, label='pump', linewidth = 1.)
    # ax4.plot(times2, pump2, color=colorPump, label='pump', linewidth = 1.)
    ax4.plot(timesNB, pump1 * 2. / 3., color=colorPump, label='pump', linewidth=1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize, labelpad=3.)
    ax2.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize, labelpad=3.)
    ax4.set_ylabel(r"$F(t) / F_0$", fontsize=fontsize, labelpad = 0.)

    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize=fontsize)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim(10., 14.)
    ax2.set_xlim(10., 14.)
    ax4.set_xlim(10., 14.)

    ax4.set_xticks([10, 11, 12, 13, 14])
    ax4.set_xticklabels(["$10$", "$11$", "$12$", "$13$", "$14$"], fontsize = fontsize)

    ax1.set_ylim(0.110, 0.142)
    ax1.set_yticks([0.12, 0.14])
    ax1.set_yticklabels(["$0.12$", "$0.14$"], fontsize = fontsize)

    ax2.set_ylim(0.110, 0.142)
    ax2.set_yticks([0.12, 0.14])
    ax2.set_yticklabels(["$0.12$", "$0.14$"], fontsize = fontsize)
    #
    ax4.set_ylim(-0.175, 0.175)
    ax4.set_yticks([-0.1, 0., 0.1])
    ax4.set_yticklabels(["$-0.1$", "$0$", "$0.1$"], fontsize = fontsize)

    legend1 = ax1.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(.0, 1.), edgecolor='black', ncol=3)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend2 = ax2.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(.0, 1.), edgecolor='black', ncol=3)
    legend2.get_frame().set_alpha(0.)
    legend2.get_frame().set_boxstyle('Square', pad=0.1)
    legend2.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(0.6, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax1.text(10., 0.1445, r"$\rm{\textbf{b.)}}$", fontsize=10, color='black', alpha=1., bbox=boxProps)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    #plt.show()

    plt.savefig('functionOfWPQuadConvNtNB.png', format='png', bbox_inches='tight', dpi=600)


def plotConvInDOccLin():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileNB1 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD300NB4TS40WPh1230.hdf5", 'r')
    fileNB2 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD300NB6TS40WPh1230.hdf5", 'r')
    fileNB3 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD300NB8TS40WPh1230.hdf5", 'r')
    fileNB4 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD300NB10TS40WPh1230.hdf5", 'r')
    fileNB5 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD300NB12TS40WPh1230.hdf5", 'r')

    #fileNB1 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB4.hdf5", 'r')
    #fileNB2 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB6TS80.hdf5", 'r')
    #fileNB3 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB8TS80.hdf5", 'r')
    #fileNB4 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB10TS80.hdf5", 'r')
    #fileNB5 = h5py.File("../data/tEvol2PhGPH50WP50WD200FD200NB12TS80.hdf5", 'r')

    # read in stuff
    fileNT1 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD300NB8TS10WPh1230.hdf5", 'r')
    fileNT2 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD300NB8TS20WPh1230.hdf5", 'r')
    fileNT3 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD300NB8TS40WPh1230.hdf5", 'r')
    fileNT4 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD300NB8TS80WPh1230.hdf5", 'r')
    fileNT5 = h5py.File("../data/tEvol2PhGPH50WP14WD200FD300NB8TS160WPh1230.hdf5", 'r')

    # readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP4 = (fileNB3['wPh'][()])[0]
    wPtWP4 = (fileNB3['wPt'][()])[0]
    wPWP4 = (fileNB3['wP'][()])[0]

    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP4, wPtWP4, wPWP4)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP4, wPtWP4, wPWP4)

    wDrive = (fileNB1['wDrive'][()])[0]

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))
    print('')
    print("W-Drive = {}".format(wDrive))

    timesNB = fileNB1['times'][()]
    timesNB = timesNB / (2. * np.pi / wDrive)
    pump1 = fileNB1['pump'][()]

    dOccNB1 = fileNB1['dOcc'][()]
    dOccNB2 = fileNB2['dOcc'][()]
    dOccNB3 = fileNB3['dOcc'][()]
    dOccNB4 = fileNB4['dOcc'][()]
    dOccNB5 = fileNB5['dOcc'][()]

    times1NT = fileNT1['times'][()]
    times1NT = times1NT / (2. * np.pi / wDrive)
    times2NT = fileNT2['times'][()]
    times2NT = times2NT / (2. * np.pi / wDrive)
    times3NT = fileNT3['times'][()]
    times3NT = times3NT / (2. * np.pi / wDrive)
    times4NT = fileNT4['times'][()]
    times4NT = times4NT / (2. * np.pi / wDrive)
    times5NT = fileNT5['times'][()]
    times5NT = times5NT / (2. * np.pi / wDrive)

    dOccNT1 = fileNT1['dOcc'][()]
    dOccNT2 = fileNT2['dOcc'][()]
    dOccNT3 = fileNT3['dOcc'][()]
    dOccNT4 = fileNT4['dOcc'][()]
    dOccNT5 = fileNT5['dOcc'][()]

    fig = plt.figure()
    fig.set_size_inches(3., 3.)
    # fig.set_size_inches(5., 5.)

    gs = gridspec.GridSpec(3, 1, height_ratios=[2, 2, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax4 = fig.add_subplot(gs[2], sharex=ax1)

    ax1.tick_params(direction='inout', length=4, width=.8)
    ax2.tick_params(direction='inout', length=4, width=.8)
    ax4.tick_params(direction='inout', length=4, width=.8)

    linewidth = 1.2

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax1.plot(timesNB, dOccNB1 + 0.5, color=color1, label=r'$N_B = 4$', linewidth=linewidth, zorder=2)
    ax1.plot(timesNB, dOccNB2 + 0.5, color=color2, label=r'$N_B = 6$', linewidth=linewidth, zorder=3)
    ax1.plot(timesNB, dOccNB3 + 0.5, color=color3, label=r'$N_B = 8$', linewidth=linewidth, zorder=4)
    ax1.plot(timesNB, dOccNB4 + 0.5, color=color4, label=r'$N_B = 10$', linewidth=linewidth, zorder=5)
    ax1.plot(timesNB, dOccNB5 + 0.5, color=color5, label=r'$N_B = 12$', linewidth=linewidth, zorder=5, linestyle='--',
             dashes=[4, 4])

    ax2.plot(times1NT, dOccNT1 + 0.5, color=color1, label=r'$N_t = 10$', linewidth=linewidth)
    ax2.plot(times2NT, dOccNT2 + 0.5, color=color2, label=r'$N_t = 20$', linewidth=linewidth)
    ax2.plot(times3NT, dOccNT3 + 0.5, color=color3, label=r'$N_t = 40$', linewidth=linewidth)
    ax2.plot(times4NT, dOccNT4 + 0.5, color=color4, label=r'$N_t = 80$', linewidth=linewidth)
    ax2.plot(times5NT, dOccNT5 + 0.5, color=color5, label=r'$N_t = 160$', linewidth=linewidth, linestyle='--',
             dashes=[4, 4])

    # ax4.plot(times1, pump1, color=colorPump, label='pump', linewidth = 1.)
    # ax4.plot(times2, pump2, color=colorPump, label='pump', linewidth = 1.)
    ax4.plot(timesNB, pump1 * 2. / 3., color=colorPump, label='pump', linewidth=1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize, labelpad=3.)
    ax2.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize, labelpad=3.)
    ax4.set_ylabel(r"$F(t) / F_0$", fontsize=fontsize, labelpad = 0.)

    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize=fontsize)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim(10., 14.)
    ax2.set_xlim(10., 14.)
    ax4.set_xlim(10., 14.)

    ax4.set_xticks([10, 11, 12, 13, 14])
    ax4.set_xticklabels(["$10$", "$11$", "$12$", "$13$", "$14$"], fontsize = fontsize)


    ax1.set_ylim(0.106, 0.131)
    ax1.set_yticks([0.11, 0.13])
    ax1.set_yticklabels(["$0.11$", "$0.13$"], fontsize = fontsize)

    ax2.set_ylim(0.106, 0.131)
    ax2.set_yticks([0.11, 0.13])
    ax2.set_yticklabels(["$0.11$", "$0.13$"], fontsize = fontsize)

    ax4.set_ylim(-0.175, 0.175)
    ax4.set_yticks([-0.1, 0., 0.1])
    ax4.set_yticklabels(["$-0.1$", "$0$", "$0.1$"], fontsize = fontsize)

    legend1 = ax1.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(.0, 1.), edgecolor='black', ncol=3)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend2 = ax2.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(.0, 1.), edgecolor='black', ncol=3)
    legend2.get_frame().set_alpha(0.)
    legend2.get_frame().set_boxstyle('Square', pad=0.1)
    legend2.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(0.6, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)


    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax1.text(10., 0.133, r"$\rm{\textbf{a.)}}$", fontsize=10, color='black', alpha=1., bbox=boxProps)


    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    #plt.show()

    plt.savefig('functionOfWPLinConvNtNB.png', format='png', bbox_inches='tight', dpi=600)


