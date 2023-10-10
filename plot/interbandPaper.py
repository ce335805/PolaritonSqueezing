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
    filenameGS = "./../data/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE10.hdf5"

    gE = 0.1

    ### get GS data ###
    wPhArr, dOcc0, dOcc1, dOccUpDn, dOccSigSig, n0, n1 = getDataFromFile(filenameGS)
    dOccIntra = dOcc0 - 2. * n0 * n0# + dOcc1 - 2. * n1 * n1
    dOccInter = dOccUpDn - 2 * n0 * n1 + dOccSigSig - 2. * n0 * n1

    print(wPhArr.shape)


    nWD = 1000
    nPhMax = np.zeros(nWD)
    nPhMean = np.zeros(nWD)
    dOccIntraMax = np.zeros(nWD)
    dOccIntraMean = np.zeros(nWD)
    #dOccInterMax = np.zeros(nWD)
    wDArr = np.zeros(nWD)
    for wDInd, wD in enumerate(wDArr):
        wDArr[wDInd] = (wDInd) / 100. + 5.
    ### get driving data ###
    for wDInd, wD in enumerate(wDArr):
        #filename = "./../data/clusterData/data/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB8TS20WPh{}.hdf5".format(int(np.rint(wD * 100)), int(np.rint(wD * 1000)))
        #filename = "./../data/gsProp2BandstHop10Ua0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wD * 100)), int(np.rint(wD * 1000)))
        filename = "./../data/clusterData/nb2Data/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wD * 100)), int(np.rint(wD * 1000)))
        #filename = "./../data/clusterData/data/gsProp2Bands1BosUa0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wD * 100)), int(np.rint(wD * 1000)))
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
        nPh = file['Nph2'][()]

        lowerTimeBound = 300
        upperTimeBound = 1300

        #nPhMax[wDInd] = np.amax(nPh)
        nPhMean[wDInd] = np.mean(nPh[lowerTimeBound : upperTimeBound])
        dOccIntraMax[wDInd] = np.amax(dOcc0 - 2. * n0 * n0)# + dOcc1 - 2. * n1 * n1)
        print(dOcc0.shape)

        dOccIntraMean[wDInd] = np.mean(dOcc0[lowerTimeBound:upperTimeBound] - 2. * n0[lowerTimeBound:upperTimeBound] * n0[lowerTimeBound:upperTimeBound])
        #dOccInterMax[wDInd] = np.amin(dOccUpDn - 2. * n0 * n1 + dOccSigSig - 2. * n0 * n1)

    print(dOccIntraMean.shape)
    print("wD Driven array shape: {}".format(wDArr.shape))


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    fig.set_size_inches(2.5, 1.6)

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    ax.plot(wPhArr, dOccIntra, color="peru", linewidth=1., label=r"GS")
    #ax.plot(wPhArr, 10000 * dOccIntra, color="peru", alpha = 0.4, linewidth=1., label=r"GS $\times 10^{4}$")
    ax.plot(wDArr, dOccIntraMean, color="teal", linewidth=1., label=r"Driven")
    #ax.plot(wPhArr, dOccIntraMax, color=cmapPink(0.3), linewidth=1., label=r"Driving - Max")

    ax2.plot(wDArr, nPhMean, color=cmapPink(0.7), linewidth=1., label=r"$\overline{N}_{\rm Bos}$")
    #ax2.plot(wPhArr, nPhMax, color=cmapPink(0.7), linewidth=1., label=r"$\overline{N_{\rm bos}}$")

    ax.axvline(10., color = 'gray', lw = 0.5)

    ax.set_yscale('log')

    ax.set_ylabel(r"$\rm{Doublon \,\, Correlations}$" + r"$\, \,\left(C \right)$", fontsize = fontsize)
    ax2.set_ylabel(r"$N_{\rm Bos}$", fontsize = fontsize)
    ax.set_xlabel(r"$\Omega \, [\Delta E]$", fontsize = fontsize)

    #ax.text(6, 0.004, r"$\times 10^{4}$", fontsize = 8, color = 'peru', alpha = 0.4)

    ax2.set_ylim(- 0.55 / 100., 0.55)
    #ax.set_ylim(- 0.026 / 100., 0.026)
    ax.set_xlim(5.2, 15.)

    ax.set_xticks([6, 10, 14])
    ax.set_xticklabels([r"$0.6$", r"$1$", "$1.4$"])

    #ax.set_yticks([0., 0.01, 0.02])
    #ax.set_yticklabels([r"$0$", r"$0.01$", r"$0.02$"])

    ax2.set_yticks([0., 0.2, 0.4])
    ax2.set_yticklabels([r"$0$", r"$0.2$", r"$0.4$"])


    legend = ax.legend(fontsize=fontsize, loc='upper left', bbox_to_anchor=(0.0, 0.55), edgecolor='black', ncol=1, handlelength = 1.5)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax2.legend(fontsize=fontsize, loc='upper right', bbox_to_anchor=(1.0, 0.9), edgecolor='black', ncol=1, handlelength = 1.5)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)
        ax2.spines[axis].set_linewidth(0.5)

    ax.text(-0.25, 1.07, r"$\rm{b.)}$", transform = ax.transAxes, fontsize = 10)

    plt.savefig('./savedPlots/timeAv.png', format='png', bbox_inches='tight', dpi=600)


def plotTimeEvolCorr():

    wArr = [9.5]

    #nT = 400
    nT = 4800

    Corr0Arr = np.zeros((len(wArr), nT))
    pumpArr = np.zeros((len(wArr), nT))
    tArr = np.zeros((len(wArr), nT))
    nBosArr = np.zeros((len(wArr), nT))
    avArr = np.zeros((len(wArr), nT))

    for wInd, wVal in enumerate(wArr):

        #filename = "./../data/gsProp2BandstHop10Ua0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wVal * 100)), int(np.rint(wVal * 1000)))
        #filename = "./../data/clusterData/nbData/dataDownload/gsProp2BandstHop10Ua0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wVal * 100)), int(np.rint(wVal * 1000)))
        filename = "./../data/clusterData/nb2Data/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wVal * 100)), int(np.rint(wVal * 1000)))
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

        tArr[wInd, :] = times
        pumpArr[wInd, :] = pump
        Corr0Arr[wInd, :] = dOcc0 - 2. * n0 * n0
        #Corr0Arr[wInd, :] = n0
        nBosArr[wInd, :] = nPh

    #compute averaged correlations
    for avInd, avVal in enumerate(avArr[0, :]):
        avArr[0, avInd] = np.sum(Corr0Arr[0, 0 : avInd]) / (avInd + 1.)


    nrow = 2
    ncol = 1
    fig = plt.figure(figsize=(7.04, 1.1), dpi=800)

    gs = gridspec.GridSpec(nrow, ncol, height_ratios=[5, 2], hspace=0.15, top=0.9, bottom=0., left=0.2, right=0.96)

    ax1 = plt.subplot(gs[0, 0])
    ax1NBos = ax1.twinx()
    ax1Drive = plt.subplot(gs[1, 0])

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    ind1 = 0
    print(wArr[ind1])
    ax1NBos.plot(tArr[ind1, :] * wArr[ind1] / 2. / np.pi, nBosArr[ind1, :], color=cmapBone(0.3), linewidth=1., label=r"$N_{\rm bos}$")
    ax1Drive.plot(tArr[ind1, :] * wArr[ind1] / 2. / np.pi, pumpArr[ind1, :], color=cmapPink(0.3), linewidth=1., label=r"$\mathrm{pump}$")
    ax1.plot(tArr[ind1, :] * wArr[ind1] / 2. / np.pi, Corr0Arr[ind1, :], color=cmapPink(0.7), linewidth=1., label=r"$C$")
    ax1.plot(tArr[ind1, :] * wArr[ind1] / 2. / np.pi, avArr[ind1, :], color='red', linewidth=1., label=r"$ \frac{1}{t} \int_{0}^{t} C(t') \, \mathrm{d}t'$")

    ax1.axhline(0., color = 'black', lw = 0.5)
    ax1Drive.axvline(tArr[0,300] * wArr[0] / 2. / np.pi, color = 'black', lw = 0.8)
    ax1Drive.axvline(tArr[0,1300] * wArr[0] / 2. / np.pi, color = 'black', lw = 0.8)

    print("t_0 in driving periods = {}".format(tArr[0,300] * wArr[0] / 2. / np.pi))
    print("t_1 in driving periods = {}".format(tArr[0,1300] * wArr[0] / 2. / np.pi))

    ax1.set_ylabel(r"$C$")
    ax1NBos.set_ylabel(r"$N_{\rm bos}$")
    ax1Drive.set_ylabel(r"$F(t) \sin(\omega_{\rm D}t)$", labelpad=8)

    ax1Drive.set_ylim(-2.5, 2.5)
    ax1Drive.set_yticks([-1, 1])

    ax1.set_xlim(tArr[ind1, 0] * wArr[ind1] / 2. / np.pi, tArr[ind1, -1] * wArr[ind1] / 2. / np.pi,)
    ax1Drive.set_xlim(tArr[ind1, 0] * wArr[ind1] / 2. / np.pi, tArr[ind1, -1] * wArr[ind1] / 2. / np.pi, )

    ax1.set_ylim(-0.0, 0.06)
    ax1NBos.set_ylim(0.0, 0.35)
    ax1Drive.set_ylim(-2.1, 2.1)

    ax1.set_xticks([])
    ax1Drive.set_xlabel(r"$t[\frac{2 \pi}{\omega_{\rm D}}]$")

    ax1Drive.set_xticks([0, 15, 50, 65, 100, 150, 200])
    ax1Drive.set_xticklabels([r"$0$", r"$t_0$", r"$50$", r"$t_1$", r"$100$", r"$150$", r"$200$"])

    legend = ax1.legend(fontsize=fontsize + 2, loc='upper left', bbox_to_anchor=(.05, 1.5), edgecolor='black', ncol=2)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax1NBos.legend(fontsize=fontsize + 2, loc='upper right', bbox_to_anchor=(1.0, 1.45), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax1Drive.legend(fontsize=fontsize + 1, loc='upper right', bbox_to_anchor=(1., 1.4), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    ax1.text(0., 1.07, r"$\rm{b.)}$", transform = ax1.transAxes, fontsize = 10)
    ax1.text(0.6, 1.07, r"$\omega_{\mathrm{D}} = 0.95 \Delta E$", transform = ax1.transAxes, fontsize = 10)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(0.5)
        ax1NBos.spines[axis].set_linewidth(0.5)
        ax1Drive.spines[axis].set_linewidth(0.5)

    plt.savefig('./savedPlots/timeEvolCorrW095.png', format='png', bbox_inches='tight', dpi=600)




def plotTimeEvolCorrSideBySide():

    wArr = [9., 10.]

    #nT = 400
    nT = 4800

    Corr0Arr = np.zeros((len(wArr), nT))
    pumpArr = np.zeros((len(wArr), nT))
    tArr = np.zeros((len(wArr), nT))
    nBosArr = np.zeros((len(wArr), nT))
    avArr = np.zeros((len(wArr), nT))

    for wInd, wVal in enumerate(wArr):

        #filename = "./../data/gsProp2BandstHop10Ua0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wVal * 100)), int(np.rint(wVal * 1000)))
        #filename = "./../data/clusterData/nbData/dataDownload/gsProp2BandstHop10Ua0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wVal * 100)), int(np.rint(wVal * 1000)))
        filename = "./../data/clusterData/nb2Data/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE1WD{}FD500NB6TS20WPh{}.hdf5".format(int(np.rint(wVal * 100)), int(np.rint(wVal * 1000)))
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

        tArr[wInd, :] = times
        pumpArr[wInd, :] = pump
        Corr0Arr[wInd, :] = dOcc0 - 2. * n0 * n0
        #Corr0Arr[wInd, :] = n0
        nBosArr[wInd, :] = nPh

    #compute averaged correlations
    for avInd, avVal in enumerate(avArr[0, :]):
        avArr[0, avInd] = np.sum(Corr0Arr[0, 0 : avInd]) / (avInd + 1.)
    for avInd, avVal in enumerate(avArr[1, :]):
        avArr[1, avInd] = np.sum(Corr0Arr[1, 0 : avInd]) / (avInd + 1.)



    nrow = 2
    ncol = 2
    fig = plt.figure(figsize=(7.04, 1.5), dpi=800)

    gs = gridspec.GridSpec(nrow, ncol, width_ratios=[1, 1], height_ratios=[5, 2],
                           hspace=0.15, wspace=0.05, top=0.9, bottom=0., left=0.2, right=0.96)

    ax1 = plt.subplot(gs[0, 0])
    ax1NBos = ax1.twinx()
    ax1Drive = plt.subplot(gs[1, 0])
    ax2 = plt.subplot(gs[0, 1])
    ax2NBos = ax2.twinx()
    ax2Drive = plt.subplot(gs[1, 1])


    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    ind1 = 0
    print(wArr[ind1])
    ax1NBos.plot(tArr[ind1, :] * wArr[ind1] / 2. / np.pi, nBosArr[ind1, :], color=cmapBone(0.3), linewidth=1., label=r"$N_{\rm bos}$")
    ax1Drive.plot(tArr[ind1, :] * wArr[ind1] / 2. / np.pi, pumpArr[ind1, :], color=cmapPink(0.3), linewidth=0.5, label=r"$\mathrm{pump}$")
    ax1.plot(tArr[ind1, :] * wArr[ind1] / 2. / np.pi, Corr0Arr[ind1, :], color=cmapPink(0.7), linewidth=1., label=r"$C$")
    ax1.plot(tArr[ind1, :] * wArr[ind1] / 2. / np.pi, avArr[ind1, :], color='red', linewidth=1., label=r"$ \frac{1}{t} \int_{0}^{t} C(t') \, \mathrm{d}t'$")

    ind2 = 1
    print(wArr[ind2])
    ax2NBos.plot(tArr[ind2, :] * wArr[ind2] / 2. / np.pi, nBosArr[ind2, :], color=cmapBone(0.3), linewidth=1., label=r"$N_{\rm bos}$")
    ax2Drive.plot(tArr[ind2, :] * wArr[ind2] / 2. / np.pi, pumpArr[ind2, :], color=cmapPink(0.3), linewidth=0.5, label=r"$\mathrm{pump}$")
    ax2.plot(tArr[ind2, :] * wArr[ind2] / 2. / np.pi, Corr0Arr[ind2, :], color=cmapPink(0.7), linewidth=1., label=r"$C$")
    ax2.plot(tArr[ind2, :] * wArr[ind2] / 2. / np.pi, avArr[ind2, :], color='red', linewidth=1., label=r"$ \frac{1}{t} \int_{0}^{t} C(t') \, \mathrm{d}t'$")


    ax1Drive.axvline(tArr[0,300] * wArr[0] / 2. / np.pi, color = 'black', lw = 0.8)
    ax1Drive.axvline(tArr[0,1300] * wArr[0] / 2. / np.pi, color = 'black', lw = 0.8)

    ax2Drive.axvline(tArr[0,300] * wArr[0] / 2. / np.pi, color = 'black', lw = 0.8)
    ax2Drive.axvline(tArr[0,1300] * wArr[0] / 2. / np.pi, color = 'black', lw = 0.8)

    print("t_0 in driving periods = {}".format(tArr[0,300] * wArr[0] / 2. / np.pi))
    print("t_1 in driving periods = {}".format(tArr[0,1300] * wArr[0] / 2. / np.pi))

    ax1.set_ylabel(r"$C$")
    ax1Drive.set_ylabel(r"$F(t) \sin(\omega_{\rm D}t)$", labelpad=8)
    ax2NBos.set_ylabel(r"$N_{\rm bos}$")


    ax1Drive.set_ylim(-2.5, 2.5)
    ax1Drive.set_yticks([-1, 1])
    ax2Drive.set_ylim(-2.5, 2.5)
    ax2Drive.set_yticks([-1, 1])

    ax1.set_xlim(tArr[ind1, 0] * wArr[ind1] / 2. / np.pi, tArr[ind1, -1] * wArr[ind1] / 2. / np.pi,)
    ax1Drive.set_xlim(tArr[ind1, 0] * wArr[ind1] / 2. / np.pi, tArr[ind1, -1] * wArr[ind1] / 2. / np.pi, )
    ax2.set_xlim(tArr[ind1, 0] * wArr[ind1] / 2. / np.pi, tArr[ind1, -1] * wArr[ind1] / 2. / np.pi,)
    ax2Drive.set_xlim(tArr[ind1, 0] * wArr[ind1] / 2. / np.pi, tArr[ind1, -1] * wArr[ind1] / 2. / np.pi, )

    ax1.set_ylim(-0.0, 0.06)
    ax1NBos.set_ylim(0.0, 0.35)
    ax1Drive.set_ylim(-2.1, 2.1)
    ax2.set_ylim(-0.0, 0.06)
    ax2NBos.set_ylim(0.0, 0.35)
    ax2Drive.set_ylim(-2.1, 2.1)

    ax1.set_xticks([])
    ax2.set_xticks([])
    ax1Drive.set_xlabel(r"$t[\frac{2 \pi}{\omega_{\rm D}}]$")
    ax2Drive.set_xlabel(r"$t[\frac{2 \pi}{\omega_{\rm D}}]$")

    ax1Drive.set_xticks([0, 15, 50, 65, 100, 150, 200])
    ax1Drive.set_xticklabels([r"$0$", r"$t_0$", r"$50$", r"$t_1$", r"$100$", r"$150$", r"$200$"])

    ax2Drive.set_xticks([0, 15, 50, 65, 100, 150, 200])
    ax2Drive.set_xticklabels([r"$0$", r"$t_0$", r"$50$", r"$t_1$", r"$100$", r"$150$", r"$200$"])

    ax2.set_yticks([])
    ax1NBos.set_yticks([])
    ax2Drive.set_yticks([])

    legend = ax1.legend(fontsize=fontsize + 2, loc='upper left', bbox_to_anchor=(.1, 1.4), edgecolor='black', ncol=2)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax2NBos.legend(fontsize=fontsize + 2, loc='upper left', bbox_to_anchor=(.0, 1.35), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend = ax1Drive.legend(fontsize=fontsize + 1, loc='upper right', bbox_to_anchor=(1., 1.3), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    #ax1.text(0., 1.07, r"$\rm{b.)}$", transform = ax1.transAxes, fontsize = 10)
    ax1.text(0.6, 0.75, r"$\omega_{\mathrm{D}} = 0.9\Delta E$", transform = ax1.transAxes, fontsize = 10)
    ax2.text(0.6, 0.75, r"$\omega_{\mathrm{D}} = \Delta E$", transform = ax2.transAxes, fontsize = 10)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(0.5)
        ax1NBos.spines[axis].set_linewidth(0.5)
        ax1Drive.spines[axis].set_linewidth(0.5)
        ax2.spines[axis].set_linewidth(0.5)
        ax2NBos.spines[axis].set_linewidth(0.5)
        ax2Drive.spines[axis].set_linewidth(0.5)

    plt.savefig('./savedPlots/timeEvolCorrSide.png', format='png', bbox_inches='tight', dpi=600)


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

    #ax1.text(400, 0.01, r"$\Omega = {}$".format(wArr[ind1]), fontsize = 8, color = 'black')
    #ax2.text(400, 0.01, r"$\Omega = {}$".format(wArr[ind2]), fontsize = 8, color = 'black')



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

