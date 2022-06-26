import numpy as np
import matplotlib.pyplot as plt
import h5py

import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors
import h5py

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

    return (wPh, dOcc0, dOcc1, dOccUpDn, dOccSigSig)



def dOccAsOfWp(fileNames):

    print("plotting Tc as function of Q")

    gEArr = np.array([0., 0.1, 0.5, 1., 2.])

    wPhs, _, _, _, _ = getDataFromFile(fileNames[0])
    nW = len(wPhs)
    print("nW = {}".format(nW))

    dOcc0Arrs = np.zeros((len(fileNames), nW))
    dOcc1Arrs = np.zeros((len(fileNames), nW))
    dOccUpDnArrs = np.zeros((len(fileNames), nW))
    dOccSigSigArrs = np.zeros((len(fileNames), nW))

    for fileInd, fileName in enumerate(fileNames):
        _, dOcc0Arrs[fileInd, :], dOcc1Arrs[fileInd, :], dOccUpDnArrs[fileInd, :], dOccSigSigArrs[fileInd, :] = getDataFromFile(fileName)

    print("dOcc0.shape = {}".format(dOcc0Arrs.shape))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3., 2.)

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    for gEInd, gE in enumerate(gEArr):

        color1 = cmapBone(gEInd / (len(gEArr) + 1) + 0.1)
        color2 = cmapPink(gEInd / (len(gEArr) + 1) + 0.1)

        ax.plot(wPhs, dOcc0Arrs[gEInd, :], color = color1)
        ax.plot(wPhs, dOcc0Arrs[gEInd, :] + dOcc1Arrs[gEInd, :], color = color2)

    ax.set_ylabel(r"$dOcc_{\rm A}$")
    ax.set_xlabel(r"$\omega_{\rm ph}$")

    plt.savefig('./savedPlots/dOccAsOfWph.png', format='png', bbox_inches='tight', dpi=600)

