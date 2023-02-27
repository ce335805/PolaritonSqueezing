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

fontsize = 8

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

def makeLevelPlot(filename):
    print("plotting nice levels")

    file = h5py.File(filename, 'r')

    gArr = file['gE'][()]
    spectrum = file['spectrum'][()]
    nc0 = file['nc0'][()]
    nd0 = file['nd0'][()]
    nc1 = file['nc1'][()]
    nd1 = file['nd1'][()]
    nBos0 = file['nBos0'][()]
    nBos1 = file['nBos1'][()]
    print("spectrum.shape = {}".format(spectrum.shape))

    spectrum = np.reshape(spectrum, (len(gArr), 576))
    print("spectrum.shape = {}".format(spectrum.shape))

    nc0 = np.reshape(nc0, (len(gArr), 576))
    nd0 = np.reshape(nd0, (len(gArr), 576))
    nc1 = np.reshape(nc1, (len(gArr), 576))
    nd1 = np.reshape(nd1, (len(gArr), 576))
    nBos0 = np.reshape(nBos0, (len(gArr), 576))
    nBos1 = np.reshape(nBos1, (len(gArr), 576))


    nrow = 2
    ncol = 1
    fig = plt.figure(figsize=(1.7, 2.), dpi=800)

    gs = gridspec.GridSpec(nrow, ncol, height_ratios=[3, 1],
                           wspace=0.35, hspace=0.1, top=0.9, bottom=0.17, left=0.16, right=0.995)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    #ax2 = ax.twinx()

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    n = 0
    ax1.plot(gArr / 10., spectrum[:, 0 + 4 * n] / 10., color = cmapBone(0.2), lw = 1.)
    ax1.plot(gArr / 10., spectrum[:, 2 + 4 * n] / 10., color = cmapBone(0.6), lw = 1.)
    ax1.plot(gArr / 10., spectrum[:, 4] / 10., color = cmapPink(0), lw = 1.)
    ax1.plot(gArr / 10., spectrum[:, 8] / 10., color = cmapPink(0.4), lw = 1.)
    ax1.plot(gArr / 10., spectrum[:, 10] / 10., color = cmapPink(0.6), lw = 1.)

    n = 0
    ax2.plot(gArr / 10., spectrum[:, 0 + 4 * n] / 10., color = cmapBone(0.2), lw = 1.)
    ax2.plot(gArr / 10., spectrum[:, 2 + 4 * n] / 10., color = cmapBone(0.6), lw = 1.)
    ax2.plot(gArr / 10., spectrum[:, 4] / 10., color = cmapPink(0), lw = 1.)
    ax2.plot(gArr / 10., spectrum[:, 8] / 10., color = cmapPink(0.4), lw = 1.)
    ax2.plot(gArr / 10., spectrum[:, 10] / 10., color = cmapPink(0.6), lw = 1.)



    ax1.set_xlim(gArr[0] / 10., 0.21)
    ax2.set_xlim(gArr[0] / 10., 0.21)
    ax1.set_ylim(0.22, 0.55)
    ax2.set_ylim(-0.1, 0.01)

    d = .02
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs, lw=0.5)  # top-left diagonal
    shift = 0.06
    ax2.plot((-d, +d), (-d - shift, +d - shift), **kwargs, lw = 0.5)  # top-right diagonal
    ax1.plot((1., 1.), (0., -0.1), **kwargs, lw = 0.5)  # top-right diagonal

    ax1.set_xticks([])

    ax2.set_xlabel(r"$g[\Delta]$", fontsize=8)
    ax1.set_ylabel(r"$E[\Delta]$", fontsize=8)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(0.5)
        ax2.spines[axis].set_linewidth(0.5)

        ax1.spines['bottom'].set_linewidth(0.)
        ax2.spines['top'].set_linewidth(0.)
        #ax1.spines['top'].set_linewidth(0.)
        #ax1.spines['right'].set_linewidth(0.)
        #ax2.spines['right'].set_linewidth(0.)


    plt.savefig('savedPlots/levelsExcited.png', format='png', bbox_inches='tight', dpi = 600)


def plotEigenstateProperties(filename):
    print("plotting nice levels")

    file = h5py.File(filename, 'r')

    gArr = file['gE'][()]
    spectrum = file['spectrum'][()]
    nc0 = file['nc0'][()]
    nd0 = file['nd0'][()]
    nc1 = file['nc1'][()]
    nd1 = file['nd1'][()]
    nBos0 = file['nBos0'][()]
    nBos1 = file['nBos1'][()]
    print("spectrum.shape = {}".format(spectrum.shape))

    spectrum = np.reshape(spectrum, (len(gArr), 576))
    print("spectrum.shape = {}".format(spectrum.shape))

    nc0 = np.reshape(nc0, (len(gArr), 576))
    nd0 = np.reshape(nd0, (len(gArr), 576))
    nc1 = np.reshape(nc1, (len(gArr), 576))
    nd1 = np.reshape(nd1, (len(gArr), 576))
    nBos0 = np.reshape(nBos0, (len(gArr), 576))
    nBos1 = np.reshape(nBos1, (len(gArr), 576))


    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax2 = ax.twinx()
    fig.set_size_inches(2.7, 1.8)

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    n = 0
    #ax.plot(gArr / 10., spectrum[:, 0 + 4 * n] / 10., color = cmapBone(0), lw = 1.5)
    #ax.plot(gArr / 10., spectrum[:, 2 + 4 * n] / 10., color = cmapBone(0.4), lw = 1.5)
    #ax.plot(gArr / 10., spectrum[:, 4] / 10., color = cmapPink(0), lw = 1.5)
    #ax.plot(gArr / 10., spectrum[:, 6] / 10., color = cmapPink(0.2), lw = 1.5)
    #ax.plot(gArr / 10., spectrum[:, 8] / 10., color = cmapPink(0.4), lw = 1.5)
    #ax.plot(gArr / 10., spectrum[:, 10] / 10., color = cmapPink(0.6), lw = 1.5)


    ax.plot(gArr, nc0[:, 11], color = 'red', lw = 1.5)
    #ax.plot(gArr, nd0[:, 11], color = 'blue', lw = 1.5)
    #ax.plot(gArr, nc1[:, 11], color = 'green', lw = 1.5)
    #ax.plot(gArr, nd1[:, 11], color = 'yellow', lw = 1.5)

    ax.plot(gArr, nBos1[:, 11], color = 'red', lw = 1.5)
    #ax.plot(gArr, nBos1[:, 1], color = 'green', lw = 1.5)

    #ax.set_xlim(gArr[0] / 10., 0.2)

    ax.set_xlabel(r"$g[\Delta]$", fontsize=8)
    ax.set_ylabel(r"$E[\Delta]$", fontsize=8)



    plt.savefig('savedPlots/eigenstateProps.png', format='png', bbox_inches='tight', dpi = 600)



