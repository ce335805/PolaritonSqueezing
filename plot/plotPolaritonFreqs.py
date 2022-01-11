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


def calcWPlus(wPh, wPt, wP):
    sqrtPart = np.sqrt((wPh ** 2 + wPt ** 2 + wP ** 2) ** 2 - 4. * wPh ** 2 * wPt ** 2)
    wPlus = 0.5 * (wPh ** 2 + wPt ** 2 + wP ** 2 + sqrtPart)

    return np.sqrt(wPlus)


def calcWMinus(wPh, wPt, wP):
    sqrtPart = np.sqrt((wPh ** 2 + wPt ** 2 + wP ** 2) ** 2 - 4. * wPh ** 2 * wPt ** 2)
    wPlus = 0.5 * (wPh ** 2 + wPt ** 2 + wP ** 2 - sqrtPart)

    return np.sqrt(wPlus)


def plotPolaritonFreqs():
    wPh = 2.
    wPt = 2.

    wPArr = np.linspace(0., 10., 100, endpoint=True) / wPt

    wPlus = calcWPlus(wPh, wPt, wPArr)
    wMinus = calcWMinus(wPh, wPt, wPArr)

    wDiff = wPlus - wMinus

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)

    fig.set_size_inches(2., 2.)

    ax.tick_params(direction='inout', length=4, width=.8)

    linewidth = 1.

    ax.plot(wPArr, wPlus / 2., color = 'olive', label = r'$\omega_+$', linewidth = linewidth)
    ax.plot(wPArr, wMinus / 2., color = 'rosybrown', label = r'$\omega_-$', linewidth = linewidth)

    #ax.set_xlabel(r"$\omega_{\rm P} / \omega_{\rm phot}$", fontsize = fontsize)
    ax.set_xlabel(r"$\rm light{-}matter$ $\rm coupling$ $[\omega_{\rm P} / \omega_{\rm phot}]$", fontsize = fontsize-1)
    ax.set_ylabel(r"$\omega / \omega_{\rm phot}$", fontsize = fontsize)

    ax.set_xlim(0., 5.)
    ax.set_ylim(0., 3.)

    ax.set_xticks([0, 2, 4])
    ax.set_xticklabels(["$0$", "$2$", "$4$"], fontsize = fontsize)
    ax.set_yticks([0., 1., 2.])
    ax.set_yticklabels(["$0$", "$1$", "$2$"], fontsize = fontsize)

    arrow = patches.FancyArrowPatch((3., calcWMinus(wPh, wPt, 3.) / 2.), (3., calcWPlus(wPh, wPt, 3.) / 2. ), arrowstyle='<->', mutation_scale=10, zorder = 100, linewidth=1., color = 'black')
    ax.add_patch(arrow)
    ax.text(3.2, 2.5 / 2., r"$\omega_{\rm P}$", fontsize = fontsize)


    legend = ax.legend(fontsize = fontsize, loc = 'upper left', bbox_to_anchor=(.0, 1.), edgecolor = 'black', ncol = 1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    plt.show()

    #plt.savefig('PolaritonFreqs.png', format='png', bbox_inches='tight', dpi = 600)
