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



def plotClassicalVSQuantumQuad():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../dataRef/tEvol2PhGPH20WP0WD200FD300NB10TS80WPh2022classic.hdf5", 'r')
    fileWP2 = h5py.File("../dataRef/tEvol2PhGPH20WP14WD200FD300NB10TS80WPh2022.hdf5", 'r')

    #fileWP1 = h5py.File("../data/tEvol2PhGPH20WP0WD200FD100NB10.hdf5", 'r')
    #fileWP2 = h5py.File("../data/tEvol2PhGPH20WP20WD200FD200NB10TS80.hdf5", 'r')

    #readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP2 = (fileWP2['wPh'][()])[0]
    wPtWP2 = (fileWP2['wPt'][()])[0]
    wPWP2 = (fileWP2['wP'][()])[0]

    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP2, wPtWP2, wPWP2)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP2, wPtWP2, wPWP2)

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

    NptWP1 = fileWP1['Npt'][()]
    NptWP2 = fileWP2['Npt'][()]

    NphWP1 = fileWP1['N1ph'][()]
    NphWP2 = fileWP2['N1ph'][()]

    fig = plt.figure()
    fig.set_size_inches(3., 3.)

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex = ax1)
    ax3 = fig.add_subplot(gs[2], sharex = ax1)
    ax4 = fig.add_subplot(gs[3], sharex = ax1)

    ax1.tick_params(direction='inout', length=4, width=.8)
    ax2.tick_params(direction='inout', length=4, width=.8)
    ax3.tick_params(direction='inout', length=4, width=.8)
    ax4.tick_params(direction='inout', length=4, width=.8)

    linewidth = 0.8

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax1.plot(times, dOccWP1 + 0.5, color=color1, label=r'Classical Drive', linewidth = linewidth)
    ax1.plot(times, dOccWP2 + 0.5, color=color3, label=r'Cavity Drive', linewidth = linewidth)

    ax2.plot(times, NphWP1, color=color1, linewidth = linewidth)
    ax2.plot(times, NphWP2, color=color3, linewidth = linewidth)

    ax3.plot(times, NptWP1, color=color1, linewidth = linewidth)
    ax3.plot(times, NptWP2, color=color3, linewidth = linewidth)

    ax4.plot(times, pump * 2/3, color=colorPump, label='pump', linewidth = 1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize, labelpad=3.)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize = fontsize, labelpad=7.)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize = fontsize, labelpad=7.)
    ax4.set_ylabel(r"$F(t) / F_0$", fontsize = fontsize, labelpad=0.)


    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize = fontsize)


    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim(1., 23)
    ax2.set_xlim(1., 23)
    ax3.set_xlim(1., 23)
    ax4.set_xlim(1., 23)
    ax1.set_xticks([5, 10, 15, 20])
    ax2.set_xticks([5, 10, 15, 20])
    ax3.set_xticks([5, 10, 15, 20])
    ax4.set_xticks([5, 10, 15, 20])
    ax4.set_xticklabels(["$5$", "$10$", "$15$", "$20$"], fontsize = fontsize)


    ax1.set_ylim(0.105, 0.145)
    ax1.set_yticks([0.12, 0.14])
    ax1.set_yticklabels(["$0.12$", "$0.14$"], fontsize = fontsize)
##
    ax2.set_ylim(-0.1, 0.8)
    ax2.set_yticks([0, 0.5])
    ax2.set_yticklabels(["$0$", "$0.5$"], fontsize = fontsize)
##
    ax3.set_ylim(-0.1, 0.8)
    ax3.set_yticks([0, 0.5])
    ax3.set_yticklabels(["$0$", "$0.5$"], fontsize = fontsize)
##
    ax4.set_ylim(-0.175, 0.175)
    ax4.set_yticks([-0.1, 0., 0.1])
    ax4.set_yticklabels(["$-0.1$", "$0$", "$0.1$"], fontsize = fontsize)


    legend1 = ax1.legend(fontsize = fontsize - 4, loc = 'upper left', bbox_to_anchor=(.0, 1.), edgecolor = 'black', ncol = 1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(0.7, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax1.text(1., 0.148, r"$\rm{\textbf{b.)}}$", fontsize=10, color='black', alpha=1., bbox=boxProps)


    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    #plt.savefig('quantumVSclassical.png', format='png', bbox_inches='tight', dpi = 600)

def plotClassicalVSQuantumLin():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileWP1 = h5py.File("../dataRef/tEvol2PhGPH50WP0WD200FD300NB10TS80WPh1230classic.hdf5", 'r')
    fileWP2 = h5py.File("../dataRef/tEvol2PhGPH50WP14WD200FD300NB10TS80WPh1230.hdf5", 'r')

    #readInPrmsAndAssert(fileWP0, fileWP01)
    wPhWP2 = (fileWP2['wPh'][()])[0]
    wPtWP2 = (fileWP2['wPt'][()])[0]
    wPWP2 = (fileWP2['wP'][()])[0]

    wMinus = plotPolaritonFreqs.calcWMinus(wPhWP2, wPtWP2, wPWP2)
    wPlus = plotPolaritonFreqs.calcWPlus(wPhWP2, wPtWP2, wPWP2)

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

    NptWP1 = fileWP1['Npt'][()]
    NptWP2 = fileWP2['Npt'][()]

    NphWP1 = fileWP1['N1ph'][()]
    NphWP2 = fileWP2['N1ph'][()]

    fig = plt.figure()
    fig.set_size_inches(3., 3.)

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex = ax1)
    ax3 = fig.add_subplot(gs[2], sharex = ax1)
    ax4 = fig.add_subplot(gs[3], sharex = ax1)

    ax1.tick_params(direction='inout', length=4, width=.8)
    ax2.tick_params(direction='inout', length=4, width=.8)
    ax3.tick_params(direction='inout', length=4, width=.8)
    ax4.tick_params(direction='inout', length=4, width=.8)

    linewidth = 0.8

    color1 = '#28384D'
    color2 = '#337343'
    color3 = '#D4AE55'
    color4 = '#734D2E'
    color5 = '#FFC3B1'
    # colorPump = '#5295A4'
    colorPump = 'red'

    ax1.plot(times, dOccWP1 + 0.5, color=color1, label=r'Classical Drive', linewidth = linewidth)
    ax1.plot(times, dOccWP2 + 0.5, color=color3, label=r'Cavity Drive', linewidth = linewidth)

    ax2.plot(times, NphWP1, color=color1, linewidth = linewidth)
    ax2.plot(times, NphWP2, color=color3, linewidth = linewidth)

    ax3.plot(times, NptWP1, color=color1, linewidth = linewidth)
    ax3.plot(times, NptWP2, color=color3, linewidth = linewidth)

    ax4.plot(times, pump * 2/3, color=colorPump, label='pump', linewidth = 1.)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize, labelpad=3.)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize = fontsize, labelpad= 7.)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize = fontsize, labelpad= 7.)
    ax4.set_ylabel(r"$F(t) / F_0$", fontsize = fontsize, labelpad=0.)


    ax4.set_xlabel(r"$t \, \, [2 \pi / \omega_{\rm Drive}]$", fontsize = fontsize)


    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim(1., 23)
    ax2.set_xlim(1., 23)
    ax3.set_xlim(1., 23)
    ax4.set_xlim(1., 23)
    ax1.set_xticks([5, 10, 15, 20])
    ax2.set_xticks([5, 10, 15, 20])
    ax3.set_xticks([5, 10, 15, 20])
    ax4.set_xticks([5, 10, 15, 20])
    ax4.set_xticklabels(["$5$", "$10$", "$15$", "$20$"], fontsize = fontsize)


    ax1.set_ylim(0.105, 0.13)
    ax1.set_yticks([0.11, 0.12])
    ax1.set_yticklabels(["$0.11$", "$0.12$"], fontsize = fontsize)
#
    ax2.set_ylim(-0.1, 1.)
    ax2.set_yticks([0, 0.5])
    ax2.set_yticklabels(["$0$", "$0.5$"], fontsize = fontsize)
#
    ax3.set_ylim(-0.1, 0.8)
    ax3.set_yticks([0, 0.5])
    ax3.set_yticklabels(["$0$", "$0.5$"], fontsize = fontsize)
#
    ax4.set_ylim(-0.175, 0.175)
    ax4.set_yticks([-0.1, 0., 0.1])
    ax4.set_yticklabels(["$-0.1$", "$0$", "$0.1$"], fontsize = fontsize)


    legend1 = ax1.legend(fontsize = fontsize - 4, loc = 'upper left', bbox_to_anchor=(.0, 1.), edgecolor = 'black', ncol = 1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend4 = ax4.legend(fontsize=fontsize - 4, loc='upper left', bbox_to_anchor=(0.7, 1.1), edgecolor='black', ncol=1)
    legend4.get_frame().set_alpha(0.)
    legend4.get_frame().set_boxstyle('Square', pad=0.1)
    legend4.get_frame().set_linewidth(0.0)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax1.text(1., 0.132, r"$\rm{\textbf{a.)}}$", fontsize=10, color='black', alpha=1., bbox=boxProps)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    #plt.savefig('quantumVSclassicalLin.png', format='png', bbox_inches='tight', dpi = 600)

