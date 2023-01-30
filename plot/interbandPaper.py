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


def getDataFromFile(fileName):
    file = h5py.File(fileName, 'r')
    wPh = file['wPh'][()]
    dOcc0 = file['dOcc0'][()]
    dOcc1 = file['dOcc1'][()]
    dOccUpDn = file['dOccUpDn'][()]
    dOccSigSig = file['dOccSigSig'][()]
    n0 = file['n0'][()]
    n1 = file['n1'][()]

    return (wPh, dOcc0, dOcc1, dOccUpDn, dOccSigSig, n0, n1)


def asOfFreq():
    print("plotting Tc as function of Q")
    filenameGS = "./../data/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE10.hdf5"
    #filenameGS = "./../data/gsProp2Bands1BosUa0Ub0Uud0Uss0epsA0epsB100gE10.hdf5"

    gE = 0.1

    ### get GS data ###
    wPhArr, dOcc0, dOcc1, dOccUpDn, dOccSigSig, n0, n1 = getDataFromFile(filenameGS)
    dOccIntra = dOcc0 - 2. * n0 * n0# + dOcc1 - 2. * n1 * n1
    dOccInter = dOccUpDn - 2 * n0 * n1 + dOccSigSig - 2. * n0 * n1

    print(dOccIntra)

    print(wPhArr.shape)


    nWD = 100
    nPhMax = np.zeros(nWD)
    nPhMean = np.zeros(nWD)
    dOccIntraMax = np.zeros(nWD)
    dOccIntraMean = np.zeros(nWD)
    #dOccInterMax = np.zeros(nWD)
    wDArr = np.zeros(nWD)
    for wDInd, wD in enumerate(wDArr):
        wDArr[wDInd] = (wDInd + 1.) / 10. + 5.
    ### get driving data ###
    for wDInd, wD in enumerate(wDArr):
        #filename = "./../data/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD100NB6TS20WPh{}.hdf5".format(int(np.rint(wD * 100)), int(np.rint(wD * 1000)))
        filename = "./../data/gsProp2BandstHop10Ua0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wD * 100)), int(np.rint(wD * 1000)))
        #print(filename)
        file = h5py.File(filename,'r')
        #times = file['times'][()]
        #pump = file['pump'][()]
        dOcc0 = file['dOcc0'][()]
        dOcc1 = file['dOcc1'][()]
        n0 = file['n0'][()]
        n1 = file['n1'][()]
        #dOccUpDn = file['dOccUpDn'][()]
        #dOccSigSig = file['dOccSigSig'][()]
        nPh = file['Nph1'][()]

        nPhMax[wDInd] = np.amax(nPh)
        nPhMean[wDInd] = np.mean(nPh)
        dOccIntraMax[wDInd] = np.amax(dOcc0 - 2. * n0 * n0)# + dOcc1 - 2. * n1 * n1)
        dOccIntraMean[wDInd] = np.mean(dOcc0[100:800] - 2. * n0[100:800] * n0[100:800])
        #dOccInterMax[wDInd] = np.amin(dOccUpDn - 2. * n0 * n1 + dOccSigSig - 2. * n0 * n1)

    print(wDArr.shape)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    fig.set_size_inches(3., 2.)

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    ax.plot(wPhArr, dOccIntra, color="peru", linewidth=1., label=r"GS")
    ax.plot(wPhArr, 10000 * dOccIntra, color="peru", alpha = 0.4, linewidth=1., label=r"GS $\times 10^{4}$")
    ax.plot(wPhArr, dOccIntraMean, color="teal", linewidth=1., label=r"Driven")

    #ax.plot(wPhArr, dOccIntraMax, color=cmapPink(0.3), linewidth=1., label=r"Driving - Max")

    ax2.plot(wPhArr, nPhMean, color=cmapPink(0.7), linewidth=1., label=r"$N_{\rm bos}$")
    #ax2.plot(wPhArr, nPhMax, color=cmapPink(0.7), linewidth=1., label=r"$\overline{N_{\rm bos}}$")

    ax.axvline(10., color = 'gray', lw = 0.5)

    ax.set_ylabel(r"$\langle n_{\uparrow} n_{\downarrow} \rangle - \langle n \rangle^2$", fontsize = fontsize)
    ax2.set_ylabel(r"$N_{\rm bos}$", fontsize = fontsize)
    ax.set_xlabel(r"$\Omega \, [\Delta]$", fontsize = fontsize)

    ax.text(6, 0.004, r"$\times 10^{4}$", fontsize = 8, color = 'peru', alpha = 0.4)

    ax2.set_ylim(- 0.55 / 100., 0.55)
    ax.set_ylim(- 0.026 / 100., 0.026)
    ax.set_xlim(5.2, 15.)

    ax.set_xticks([6, 10, 14])
    ax.set_xticklabels([r"$0.6$", r"$\Delta = \Omega$", "$1.4$"])

    ax.set_yticks([0., 0.01, 0.02])
    ax.set_yticklabels([r"$0$", r"$0.01$", r"$0.02$"])

    ax2.set_yticks([0., 0.2, 0.4])
    ax2.set_yticklabels([r"$0$", r"$0.2$", r"$0.4$"])


    legend = ax.legend(fontsize=fontsize, loc='upper left', bbox_to_anchor=(0.0, 0.65), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax2.legend(fontsize=fontsize, loc='upper right', bbox_to_anchor=(1.0, 0.6), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)
        ax2.spines[axis].set_linewidth(0.5)

    plt.savefig('./savedPlots/timeAv.png', format='png', bbox_inches='tight', dpi=600)


def plotTimeEvolCorr():

    #wArr = [6., 7., 8., 9., 10.5, 11., 12., 13., 14.]
    wArr = [6., 8.]
    wArr = [9.8, 10.1]
    #wArr = [10.3, 12.]

    #nT = 400
    nT = 2400

    Corr0Arr = np.zeros((len(wArr), nT))
    pumpArr = np.zeros((len(wArr), nT))
    tArr = np.zeros((len(wArr), nT))
    nBosArr = np.zeros((len(wArr), nT))

    for wInd, wVal in enumerate(wArr):

        filename = "./../data/gsProp2BandstHop10Ua0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD200NB6TS20WPh{}.hdf5".format(int(np.rint(wArr[wInd] * 100)), int(np.rint(wArr[wInd] * 1000)))
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

        tArr[wInd, :] = times
        pumpArr[wInd, :] = pump
        Corr0Arr[wInd, :] = dOcc0 - 2. * n0 * n0
        #Corr0Arr[wInd, :] = n0
        nBosArr[wInd, :] = nPh


    nrow = 2
    ncol = 2
    fig = plt.figure(figsize=(7.04, 2.3), dpi=800)

    gs = gridspec.GridSpec(nrow, ncol, height_ratios=[2, 1], width_ratios=[1, 1],
                           wspace=0.5, hspace=0., top=0.9, bottom=0.17, left=0.2, right=0.96)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])

    ax1NBos = ax1.twinx()
    ax2NBos = ax2.twinx()

    ax1Drive = plt.subplot(gs[1, 0])
    ax2Drive = plt.subplot(gs[1, 1])

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    ind1 = 0
    print(wArr[ind1])
    ax1NBos.plot(tArr[ind1, :] * wArr[ind1], nBosArr[ind1, :], color=cmapBone(0.3), linewidth=1., label=r"$N_{\rm bos}$")
    ax1Drive.plot(tArr[ind1, :] * wArr[ind1], pumpArr[ind1, :], color=cmapPink(0.3), linewidth=1., label=r"pump")
    ax1.plot(tArr[ind1, :] * wArr[ind1], Corr0Arr[ind1, :], color=cmapPink(0.7), linewidth=1., label=r"$\langle n_0^2 \rangle - \langle n_0 \rangle^2$")

    ind2 = 1
    print(wArr[ind2])
    ax2NBos.plot(tArr[ind2, :] * wArr[ind2], nBosArr[ind2, :], color=cmapBone(0.3), linewidth=1., label=r"$N_{\rm bos}$")
    ax2Drive.plot(tArr[ind2, :] * wArr[ind2], pumpArr[ind2, :], color=cmapPink(0.3), linewidth=1., label=r"pump")
    ax2.plot(tArr[ind2, :] * wArr[ind2], Corr0Arr[ind2, :], color=cmapPink(0.7), linewidth=1., label=r"$\langle n_0^2 \rangle - \langle n_0 \rangle^2$")


    ax1.set_ylabel(r"$\langle n_{\uparrow} n_{\downarrow} \rangle - \langle n \rangle^2$")
    ax1NBos.set_ylabel(r"$N_{\rm bos}$")
    ax1Drive.set_ylabel(r"$F(t) \sin(\omega_{\rm D}t)$")
    ax2.set_ylabel(r"$\langle n_{\uparrow} n_{\downarrow} \rangle - \langle n \rangle^2$")
    ax2NBos.set_ylabel(r"$N_{\rm bos}$")
    ax2Drive.set_ylabel(r"$F(t) \sin(\omega_{\rm D}t)$")

    #ax1.set_ylim(0., 0.15)
    #ax2.set_ylim(0., 0.15)

    #ax1NBos.set_ylim(0., 0.6)
    #ax2NBos.set_ylim(0., 0.6)

    ax1Drive.set_ylim(-2.5, 2.5)
    ax1Drive.set_yticks([-1, 0, 1])
    ax2Drive.set_ylim(-2.5, 2.5)
    ax2Drive.set_yticks([-1, 0, 1])


    ax1.set_xlim(tArr[ind1, 0] * wArr[ind1], tArr[ind1, -1] * wArr[ind1],)
    ax1Drive.set_xlim(tArr[ind1, 0] * wArr[ind1], tArr[ind1, -1] * wArr[ind1], )

    ax2.set_xlim(tArr[ind2, 0] * wArr[ind2], tArr[ind2, -1] * wArr[ind2],)
    ax2Drive.set_xlim(tArr[ind2, 0] * wArr[ind2], tArr[ind2, -1] * wArr[ind2], )

    ax1.set_xticks([])
    ax1Drive.set_xlabel(r"$t[\frac{1}{\omega_{\rm D}}]$")

    ax2.set_xticks([])
    ax2Drive.set_xlabel(r"$t[\frac{1}{\omega_{\rm D}}]$")

    ax1.text(400, 0., r"$\Omega = {}$".format(wArr[ind1]), fontsize = 8, color = 'black')
    ax2.text(400, 0., r"$\Omega = {}$".format(wArr[ind2]), fontsize = 8, color = 'black')

    legend = ax1.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(.0, 1.2), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax1NBos.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(1.0, 1.2), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax1Drive.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(1., 1.), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax2.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(.0, 1.2), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax2NBos.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(1., 1.2), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax2Drive.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(1., 1.), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(0.5)
        ax1NBos.spines[axis].set_linewidth(0.5)
        ax1Drive.spines[axis].set_linewidth(0.5)
        ax2.spines[axis].set_linewidth(0.5)
        ax2NBos.spines[axis].set_linewidth(0.5)
        ax2Drive.spines[axis].set_linewidth(0.5)

    plt.savefig('./savedPlots/timeEvolCorr3.png', format='png', bbox_inches='tight', dpi=600)

def plotTimeEvolOcc():

    #wArr = [6., 8.]
    #wArr = [9.4, 10.]
    wArr = [10.3, 12.]

    nT = 2400

    n0Arr = np.zeros((len(wArr), nT))
    n1Arr = np.zeros((len(wArr), nT))
    pumpArr = np.zeros((len(wArr), nT))
    tArr = np.zeros((len(wArr), nT))
    nBosArr = np.zeros((len(wArr), nT))

    for wInd, wVal in enumerate(wArr):

        filename = "./../data/clusterData/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE2WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wArr[wInd] * 100)), int(np.rint(wArr[wInd] * 1000)))
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

        tArr[wInd, :] = times
        pumpArr[wInd, :] = pump
        #Corr0Arr[wInd, :] = dOcc0 - 2. * n0 * n0
        n0Arr[wInd, :] = 2. * n0
        n1Arr[wInd, :] = 2. * n1
        nBosArr[wInd, :] = nPh


    nrow = 2
    ncol = 2
    fig = plt.figure(figsize=(7.04, 2.3), dpi=800)

    gs = gridspec.GridSpec(nrow, ncol, height_ratios=[2, 1], width_ratios=[1, 1],
                           wspace=0.5, hspace=0., top=0.9, bottom=0.17, left=0.2, right=0.96)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])

    ax1NBos = ax1.twinx()
    ax2NBos = ax2.twinx()

    ax1Drive = plt.subplot(gs[1, 0])
    ax2Drive = plt.subplot(gs[1, 1])

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    ind1 = 0
    print(wArr[ind1])
    ax1NBos.plot(tArr[ind1, :] * wArr[ind1], nBosArr[ind1, :], color=cmapBone(0.3), linewidth=1., label=r"$N_{\rm bos}$")
    ax1Drive.plot(tArr[ind1, :] * wArr[ind1], pumpArr[ind1, :], color=cmapPink(0.3), linewidth=1., label=r"pump")
    #ax1.plot(tArr[ind1, :] * wArr[ind1], n0Arr[ind1, :], color=cmapPink(0.7), linewidth=1., label=r"$n_0$")
    ax1.plot(tArr[ind1, :] * wArr[ind1], n1Arr[ind1, :], color=cmapPink(0.6), linewidth=1., label=r"$n_1$")

    ind2 = 1
    print(wArr[ind2])
    ax2NBos.plot(tArr[ind2, :] * wArr[ind2], nBosArr[ind2, :], color=cmapBone(0.3), linewidth=1., label=r"$N_{\rm bos}$")
    ax2Drive.plot(tArr[ind2, :] * wArr[ind2], pumpArr[ind2, :], color=cmapPink(0.3), linewidth=1., label=r"pump")
    #ax2.plot(tArr[ind2, :] * wArr[ind2], n0Arr[ind2, :], color=cmapPink(0.7), linewidth=1., label=r"$n_0$")
    ax2.plot(tArr[ind2, :] * wArr[ind2], n1Arr[ind2, :], color=cmapPink(0.6), linewidth=1., label=r"$n_1$")


    ax1.set_ylabel(r"$\langle n \rangle$")
    ax1NBos.set_ylabel(r"$N_{\rm bos}$")
    ax1Drive.set_ylabel(r"$F(t) \sin(\omega_{\rm D}t)$")
    ax2.set_ylabel(r"$\langle n \rangle$")
    ax2NBos.set_ylabel(r"$N_{\rm bos}$")
    ax2Drive.set_ylabel(r"$F(t) \sin(\omega_{\rm D}t)$")

    ax1.set_ylim(0., 0.2)
    ax2.set_ylim(0., 0.2)

    ax1NBos.set_ylim(0., 0.6)
    ax2NBos.set_ylim(0., 0.6)

    ax1Drive.set_ylim(-2.5, 2.5)
    ax1Drive.set_yticks([-1, 0, 1])
    ax2Drive.set_ylim(-2.5, 2.5)
    ax2Drive.set_yticks([-1, 0, 1])

    ax1.set_xlim(tArr[ind1, 0] * wArr[ind1], tArr[ind1, -1] * wArr[ind1],)
    ax1Drive.set_xlim(tArr[ind1, 0] * wArr[ind1], tArr[ind1, -1] * wArr[ind1], )

    ax2.set_xlim(tArr[ind2, 0] * wArr[ind2], tArr[ind2, -1] * wArr[ind2],)
    ax2Drive.set_xlim(tArr[ind2, 0] * wArr[ind2], tArr[ind2, -1] * wArr[ind2], )

    ax1.set_xticks([])
    ax1Drive.set_xlabel(r"$t[\frac{1}{\omega_{\rm D}}]$")

    ax2.set_xticks([])
    ax2Drive.set_xlabel(r"$t[\frac{1}{\omega_{\rm D}}]$")

    ax1.text(400, 0.01, r"$\Omega = {}$".format(wArr[ind1]), fontsize = 8, color = 'black')
    ax2.text(400, 0.01, r"$\Omega = {}$".format(wArr[ind2]), fontsize = 8, color = 'black')



    legend = ax1.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(.1, 1.2), edgecolor='black', ncol=2)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax1NBos.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(1., 1.2), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax1Drive.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(1., 1.), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax2.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(.1, 1.2), edgecolor='black', ncol=2)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax2NBos.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(1., 1.2), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax2Drive.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(1., 1.), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(0.5)
        ax1NBos.spines[axis].set_linewidth(0.5)
        ax1Drive.spines[axis].set_linewidth(0.5)
        ax2.spines[axis].set_linewidth(0.5)
        ax2NBos.spines[axis].set_linewidth(0.5)
        ax2Drive.spines[axis].set_linewidth(0.5)

    plt.savefig('./savedPlots/timeEvolOcc3.png', format='png', bbox_inches='tight', dpi=600)

