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
from matplotlib.patches import ConnectionPatch


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
mpl.rcParams['xtick.major.width'] = .5
mpl.rcParams['ytick.major.width'] = .5
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


def plotTimeEvolConvergenceNB():

    wVal = 10.
    nBArr = np.array([4, 6, 8, 10])


    #nT = 400
    nT = 1600

    Corr0Arr = np.zeros((len(nBArr), nT))
    pumpArr = np.zeros((len(nBArr), nT))
    tArr = np.zeros((len(nBArr), nT))
    nBosArr = np.zeros((len(nBArr), nT))

    for nBInd, nBVal in enumerate(nBArr):

        filename = "./../data/clusterData/nbData/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB{}TS20WPh{}.hdf5".format(int(wVal * 100), int(nBVal), int(wVal * 1000))
        #filename = "./../data/gsProp2BandstHop10Ua0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD200NB6TS20WPh{}.hdf5".format(int(np.rint(wArr[wInd] * 100)), int(np.rint(wArr[wInd] * 1000)))
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

        tArr[nBInd, :] = times
        pumpArr[nBInd, :] = pump
        Corr0Arr[nBInd, :] = dOcc0 - 2. * n0 * n0
        #Corr0Arr[wInd, :] = n0
        nBosArr[nBInd, :] = nPh


    nrow = 2
    ncol = 2
    fig = plt.figure(figsize=(7.04, 1.5), dpi=800)

    gs = gridspec.GridSpec(nrow, ncol, height_ratios=[2, 1], width_ratios=[1, 1],
                           wspace=0.1, hspace=0., top=0.9, bottom=0.17, left=0.2, right=0.96)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])

    ax1NBos = ax1.twinx()
    ax1Drive = plt.subplot(gs[1, 0])

    ax2NBos = ax2.twinx()
    ax2Drive = plt.subplot(gs[1, 1])


    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    for nBInd, nBVal in enumerate(nBArr):

        color1 = cmapPink((nBInd + 1.) / (len(nBArr) + 2.))
        color2 = cmapBone((nBInd + 1.) / (len(nBArr) + 2.))

        if(nBInd == len(nBArr) - 1):
            ax1NBos.plot(tArr[nBInd, :] * wVal / 2. / np.pi, nBosArr[nBInd, :], color=color2, linewidth=.7, linestyle = '--', label=r"$N_{\rm bos}, N_{\rm Max} = $" + "{}".format(nBVal))
            ax1.plot(tArr[nBInd, :] * wVal / 2. / np.pi, Corr0Arr[nBInd, :], color=color1, linewidth=.7, linestyle = '--', label=r"$C, N_{\rm Max} = $" + "{}".format(nBVal))
        else:
            ax1NBos.plot(tArr[nBInd, :] * wVal / 2. / np.pi, nBosArr[nBInd, :], color=color2, linewidth=.7, label=r"$N_{\rm bos}, N_{\rm Max} = $" + "{}".format(nBVal))
            ax1.plot(tArr[nBInd, :] * wVal / 2. / np.pi, Corr0Arr[nBInd, :], color=color1, linewidth=.7, label=r"$C, N_{\rm Max} = $" + "{}".format(nBVal))

        if(nBInd == len(nBArr) - 1):
            ax2NBos.plot(tArr[nBInd, :] * wVal / 2. / np.pi, nBosArr[nBInd, :], color=color2, linewidth=.7, linestyle = '--', label=r"$N_{\rm bos}, N_{\rm Max} = $" + "{}".format(nBVal))
            ax2.plot(tArr[nBInd, :] * wVal / 2. / np.pi, Corr0Arr[nBInd, :], color=color1, linewidth=.7, linestyle = '--', label=r"$\langle n_0^2 \rangle - \langle n_0 \rangle^2$")
        else:
            ax2NBos.plot(tArr[nBInd, :] * wVal / 2. / np.pi, nBosArr[nBInd, :], color=color2, linewidth=.7, label=r"$N_{\rm bos}, N_{\rm Max} = $" + "{}".format(nBVal))
            ax2.plot(tArr[nBInd, :] * wVal / 2. / np.pi, Corr0Arr[nBInd, :], color=color1, linewidth=.7, label=r"$\langle n_0^2 \rangle - \langle n_0 \rangle^2$")


    ax1Drive.plot(tArr[0, :] * wVal / 2. / np.pi, pumpArr[0, :], color=cmapPink(0.3), linewidth=1., label=r"pump")
    ax2Drive.plot(tArr[0, :] * wVal / 2. / np.pi, pumpArr[0, :], color=cmapPink(0.3), linewidth=1., label=r"pump")

    ax1.set_xlim(0, 80)
    ax1Drive.set_xlim(0, 80)
    ax2.set_xlim(70, 80)
    ax2Drive.set_xlim(70, 80)

    ax1.set_yticks([0, 0.025, 0.05])
    ax1.set_yticklabels([r"$0$", r"$0.025$", r"$0.05$"])

    ax1.set_ylabel(r"$C$")
    ax2NBos.set_ylabel(r"$N_{\rm bos}$")
    ax1Drive.set_ylabel(r"$F(t) \sin(\omega_{\rm D}t)$", labelpad=12)

    ax1Drive.set_xlabel(r"$t[\frac{2\pi}{\omega_{\rm D}}]$")
    ax2Drive.set_xlabel(r"$t[\frac{2\pi}{\omega_{\rm D}}]$")


    ax2.set_yticks([])
    ax1NBos.set_yticks([])
    ax2.set_yticks([])
    ax2Drive.set_yticks([])
    ax1.set_xticks([])
    ax2.set_xticks([])


    legend = ax1.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(.15, 1.4), edgecolor='black', ncol=2)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax2NBos.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(0.9, 1.4), edgecolor='black', ncol=2)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)
#
#    legend = ax1Drive.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(1., 1.), edgecolor='black', ncol=1)
#    legend.get_frame().set_alpha(0.)
#    legend.get_frame().set_boxstyle('Square', pad=0.1)
#    legend.get_frame().set_linewidth(0.0)

    ax1.text(0., 1.07, r"$\rm{a.)}$", transform = ax1.transAxes, fontsize = 10)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(0.5)
        ax1NBos.spines[axis].set_linewidth(0.5)
        ax1Drive.spines[axis].set_linewidth(0.5)
        ax2.spines[axis].set_linewidth(0.5)
        ax2NBos.spines[axis].set_linewidth(0.5)
        ax2Drive.spines[axis].set_linewidth(0.5)

    plt.savefig('./savedPlots/timeEvolConvergence.png', format='png', bbox_inches='tight', dpi=600)
