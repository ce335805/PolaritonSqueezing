import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors
import h5py

fontsize = 14

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

def avOverArr(arr, avOver):

    averagedArr = np.zeros(arr.shape)

    for ind in np.arange(len(arr)):
        for avInd in np.arange(avOver):
            averagedArr[ind] += arr[(ind - (avInd - avOver // 2)) % len(arr)] / avOver
    return averagedArr


def plotTimeEvol():
    print("plotting some beautiful time evolution")

    # read in stuff
    file = h5py.File("../data/timeEvolTwoPhonWP0N10.hdf5", 'r')
    times = file['times'][()]
    pump = file['pump'][()]
    dOcc = file['dOcc'][()]
    Xpt = file['Xpt'][()]
    XptSqr = file['XptSqr'][()]
    Npt = file['Npt'][()]
    Xph = file['X1ph'][()]
    XphSqr = file['X1phSqr'][()]
    Nph = file['N1ph'][()]

    Xph2 = file['X2ph'][()]
    XphSqr2 = file['X2phSqr'][()]
    Nph2 = file['N2ph'][()]

    print("times.shape = {}".format(times.shape))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    #ax.plot(times, dOcc + 0.5, color='rosybrown', label='dOcc, $\omega_{\rm P} = 0.2$')
    ax.plot(times, Nph, color='olive', label='$N_{ph, 1}$')
    #ax.plot(times, Nph2, color='rosybrown', label='$N_{ph, 2}$')
    #ax.plot(times, (Nph - Nph2) - (Nph[0] - Nph2[0]), color='olive', label=r'$\langle X^2 \rangle$', linewidth = 0.8)
    #ax.plot(times, (XphSqr - XphSqr2) - (XphSqr[0] - XphSqr2[0]), color='olive', label=r'$\langle X^2 \rangle$', linewidth = 0.8)
    #ax.plot(times, Xph, color='olive', label=r'$\langle X^2 \rangle$')
    #ax.plot(times, Xph2, color='rosybrown', label=r'$\langle X^2 \rangle$')

    plt.legend()

    # ax.set_ylim(-1e10, 1e10)

    plt.show()

def plotTimeEvolManyCurves():
    print("plotting some beautiful time evolution")

    # read in stuff
    fileW20 = h5py.File("../data/timeEvolOnePhonWP20N10.hdf5", 'r')
    times = fileW20['times'][()]
    pump = fileW20['pump'][()]
    dOccW20 = fileW20['dOcc'][()]
    XptW20 = fileW20['Xpt'][()]
    XptSqrW20 = fileW20['XptSqr'][()]
    NptW20 = fileW20['Npt'][()]
    XphW20 = fileW20['X1ph'][()]
    XphSqrW20 = fileW20['X1phSqr'][()]
    NphW20 = fileW20['N1ph'][()]

    fileW10 = h5py.File("../data/timeEvolOnePhonWP10N10.hdf5", 'r')
    dOccW10 = fileW10['dOcc'][()]
    XptW10 = fileW10['Xpt'][()]
    XptSqrW10 = fileW10['XptSqr'][()]
    NptW10 = fileW10['Npt'][()]
    XphW10 = fileW10['X1ph'][()]
    XphSqrW10 = fileW10['X1phSqr'][()]
    NphW10 = fileW10['N1ph'][()]

    fileW15 = h5py.File("../data/timeEvolOnePhonWP14N10.hdf5", 'r')
    dOccW15 = fileW15['dOcc'][()]
    XptW15 = fileW15['Xpt'][()]
    XptSqrW15 = fileW15['XptSqr'][()]
    NptW15 = fileW15['Npt'][()]
    XphW15 = fileW15['X1ph'][()]
    XphSqrW15 = fileW15['X1phSqr'][()]
    NphW15 = fileW15['N1ph'][()]

    fileW5 = h5py.File("../data/timeEvolOnePhonWP5N10.hdf5", 'r')
    dOccW5 = fileW5['dOcc'][()]
    XptW5 = fileW5['Xpt'][()]
    XptSqrW5 = fileW5['XptSqr'][()]
    NptW5 = fileW5['Npt'][()]
    XphW5 = fileW5['X1ph'][()]
    XphSqrW5 = fileW5['X1phSqr'][()]
    NphW5 = fileW5['N1ph'][()]

    fileW1 = h5py.File("../data/timeEvolOnePhonWP1N10.hdf5", 'r')
    dOccW1 = fileW1['dOcc'][()]
    XptW1 = fileW1['Xpt'][()]
    XptSqrW1 = fileW1['XptSqr'][()]
    NptW1 = fileW1['Npt'][()]
    XphW1 = fileW1['X1ph'][()]
    XphSqrW1 = fileW1['X1phSqr'][()]
    NphW1 = fileW1['N1ph'][()]

    print("times.shape = {}".format(times.shape))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # dOccQAv = avOverArr(dOccQ, 7)

    # ax.plot(times, pump, color = 'cornflowerblue' , label = 'pump')
    # ax.plot(times, NphC, color = 'olive' , label = r'$N_{\rm phon}$, $\omega_{\rm P} = 0.1$')
    # ax.plot(times, NphQ, color = 'rosybrown' , label = r'$N_{\rm phon}$, $\omega_{\rm P} = 0.05$')
    # ax.plot(times, NptC, color = 'mediumseagreen' , label = r'$N_{\rm phot}$')
    # ax.plot(times, 0.2 * (XphSqrC - XphSqrC[0]), label = r'$\langle X_{\rm phon}^2 \rangle$', color = 'olive')
    # ax.plot(times, XptSqrC * 0.1, label = r'$\langle X_{\rm phot}^2 \rangle$', color = 'rosybrown')
    # ax.plot(times, XphC, label = r'$\langle X_{\rm phon} \rangle$', color = 'olive')
    # ax.plot(times, XptC, label = r'$\langle X_{\rm phot} \rangle$', color = 'rosybrown')

    ax.plot(times, dOccW20 + .5, color='rosybrown', label='dOcc, $\omega_{\rm P} = 0.2$')
    ax.plot(times, dOccW15 + .5, color='c', label='dOcc, $\omega_{\rm P} = 0.15$')
    ax.plot(times, dOccW10 + .5, color='olive', label='dOcc, $\omega_{\rm P} = 0.1$')
    ax.plot(times, dOccW5 + .5, color='peru', label='dOcc, $\omega_{\rm P} = 0.05$')
    ax.plot(times, dOccW1 + .5, color='black', label='dOcc, $\omega_{\rm P} = 0.01$')
    # ax.plot(times, ((dOccC + .5) - (dOccC[0] + .5)) , color = 'rosybrown', label = 'dOcc, $\omega_{\rm P} = 0.1$')
    # ax.plot(times, ((dOccC + .5) - (dOccC[0] + .5)) / ((XphSqrC - XphSqrC[0]) + 1e-9) , color = 'rosybrown', label = 'dOcc, $\omega_{\rm P} = 0.1$')

    # ax.plot(times, (dOccC + .5) * 100 - 10.71, color = 'rosybrown', label = 'dOcc - Classical')
    # ax.plot(times, XphSqrC * 0.2, label = r'$\langle X_{\rm phon}^2 \rangle$', color = 'olive')

    plt.legend()

    # ax.set_ylim(-1e10, 1e10)

    plt.show()

