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
    filenameGS = "./../data/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE30.hdf5"

    gE = 0.3

    ### get GS data ###
    wPhArr, dOcc0, dOcc1, dOccUpDn, dOccSigSig, n0, n1 = getDataFromFile(filenameGS)
    dOccIntra = dOcc0 - 2. * n0 * n0 + dOcc1 - 2. * n1 * n1
    dOccInter = dOccUpDn - 2 * n0 * n1 + dOccSigSig - 2. * n0 * n1


    nWD = 50
    nPhMax = np.zeros(nWD)
    nPhMean = np.zeros(nWD)
    dOccIntraMax = np.zeros(nWD)
    dOccIntraMean = np.zeros(nWD)
    #dOccInterMax = np.zeros(nWD)
    wDArr = np.zeros(nWD)
    for wDInd, wD in enumerate(wDArr):
        wDArr[wDInd] = (wDInd + 1.) / 5. + 5.
    ### get driving data ###
    for wDInd, wD in enumerate(wDArr):
        filename = "./../data/gsProp2BandsUa0Ub0Uud0Uss0epsA0epsB100gE2WD{}FD100NB4TS20WPh{}.hdf5".format(int(np.rint(wD * 100)), int(np.rint(wD * 1000)))
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
        dOccIntraMax[wDInd] = np.amax(dOcc0 - 2. * n0 * n0 + dOcc1 - 2. * n1 * n1)
        dOccIntraMean[wDInd] = np.mean(dOcc0 - 2. * n0 * n0 + dOcc1 - 2. * n1 * n1)
        #dOccInterMax[wDInd] = np.amin(dOccUpDn - 2. * n0 * n1 + dOccSigSig - 2. * n0 * n1)



    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3., 2.)

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    ax.plot(wPhArr, dOccIntra, color=cmapBone(0.3), linewidth=1., label=r"$g = {}$".format(gE))
    #ax.plot(wPhArr, dOccIntraMax, color='red', linewidth=1., label=r"$g = {}$".format(gE))
    ax.plot(wPhArr, dOccIntraMean, color=cmapPink(0.3), linewidth=1., label=r"$g = {}$".format(gE))
    ax.plot(wPhArr, nPhMean, color=cmapPink(0.7), linewidth=1., label=r"$g = {}$".format(gE))

    ax.set_ylabel(r"$\langle n_{\uparrow} n_{\downarrow} \rangle - \langle n \rangle^2$")
    ax.set_xlabel(r"$\omega_{\rm ph}$")



    legend = ax.legend(fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(1.0, 1.0), edgecolor='black', ncol=1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    plt.savefig('./savedPlots/plotInterbandPaper.png', format='png', bbox_inches='tight', dpi=600)


