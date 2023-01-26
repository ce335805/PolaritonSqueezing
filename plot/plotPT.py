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



def calcPTN(nPhot, Delta, Omega):

    term1 = -(nPhot+1.)*(nPhot+2.)/(Omega * (Delta + Omega)**2)
    term2 =  nPhot * (nPhot-1.)/(Omega * (Delta - Omega)**2)
    term3 = -(nPhot**2 + 3. * nPhot + 1) / (Delta * (Delta + Omega)**2)
    term4 = - (nPhot**2 - nPhot) / (Delta * (Delta - Omega)**2)
    term5 = -2. * nPhot * (nPhot + 1.) / (Delta * (Delta**2 - Omega**2))
    term6 = - 4. * (nPhot + 1) / (Delta + Omega)**3
    term7 = 4. * nPhot / (Delta - Omega)**3

    return term1 + term2 + term3 + term4 + term5 + term6 + term7

def calcPTNM(n, m, Del, Om):
    term1 = - (n + 1.)*(n + 2.) / Om / (Del + Om)**2
    term2 = n * (n - 1.) / Om / (Del - Om)**2
    term3 = - (4. * (n + 1.)**2 - n * (n + 1.) - m * (m + 1.)) / 2. / Del / (Del + Om)**2
    term4 = - (4. * n**2 - n * (m + 1.) - m * (n + 1.)) / 2. / Del / (Del - Om)**2
    term5 = - (4. * n * (n + 1.) - n * (m + 1.) - m * (n + 1.)) / 2. / Del / (Del**2 - Om**2)
    term6 = - 4. * ((n + 1.) * (n + 2.) - (n + 1.) * (m + 1.)) / (Del + Om)**3
    term7 = - 4. * (n * (n - 1.) - n * m) / (Del - Om)**3

    terms = np.array([term1, term2, term3, term4, term5, term6, term7])
    return np.sum(terms, axis=0)


def calPTForAlpha(alpha, Del, Om):

    expectation = 0.
    for n in np.arange(10):
        for m in np.arange(10):
            expectation += alpha**(2. * n) * alpha**(2. * m) / np.math.factorial(n) / np.math.factorial(m) * calcPTNM(n, m, Del, Om)

    return expectation * np.exp(- 2. * np.abs(alpha)**2)

def plotPTN():
    Delta = 1.
    nPhot = 0.
    OmArr = np.linspace(0.1, 2.0, 100)
    ptArrN0 = calcPTN(nPhot, Delta, OmArr)

    nPhot = 1.
    ptArrN1 = calcPTN(nPhot, Delta, OmArr)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3., 2.)

    ax.plot(OmArr, ptArrN0, lw = 1., color = 'peru')
    ax.plot(OmArr, ptArrN1, lw = 1., color = 'teal')
    #ax.plot(DeltaArr, 0.95 * ptArrN0 + 0.05 * ptArrN1, lw = 1., color = 'indianred')
    ax.axhline(0., color = 'gray', lw = 0.8)

    ax.set_ylim(-30., 30.)
    ax.set_xlim(0.1, 2.)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)

    plt.savefig('./savedPlots/ptPlot.png', format='png', bbox_inches='tight', dpi=600)


def plotPTAlpha():
    Delta = 1.
    alpha = 0.
    OmArr = np.linspace(0.1, 2.0, 100)
    ptArrN0 = calPTForAlpha(alpha, Delta, OmArr)

    alpha = np.sqrt(0.3)
    print(alpha)
    ptArrAlpha = calPTForAlpha(alpha, Delta, OmArr)
    print(ptArrAlpha.shape)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3., 2.)

    ax.plot(OmArr, ptArrN0, lw = 1., color = 'peru', label = r"$|\alpha|^2 = 0$")
    ax.plot(OmArr, ptArrAlpha, lw = 1., color = 'teal', label = r"$|\alpha|^2 = 0.3$")
    #ax.plot(DeltaArr, 0.95 * ptArrN0 + 0.05 * ptArrN1, lw = 1., color = 'indianred')
    ax.axhline(0., color = 'gray', lw = 0.5)

    ax.set_ylim(- 30., 5.)
    ax.set_xlim(0.1, 2.)



    ax.set_xlabel(r"$\Omega \, [\Delta]$", fontsize = fontsize)
    ax.set_ylabel(r"$\left(E^{(4)}_{\uparrow \downarrow, -} - E^{(4)}_{\uparrow, \downarrow} \right) / g^4 \, [\Delta]$", fontsize = fontsize)

    ax.set_yticks([0., -20.])
    ax.set_yticklabels([r"$0$", r"$-20$"])

    ax.set_xticks([0.5, 1., 1.5])
    ax.set_xticklabels([r"$0.5$", r"$\Delta = \Omega$", r"$1.5$"])

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)

    legend1 = ax.legend(fontsize = fontsize - 2, loc = 'lower right', bbox_to_anchor=(1.0, 0.), edgecolor = 'black', ncol = 1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    plt.savefig('./savedPlots/ptPlotAlpha.png', format='png', bbox_inches='tight', dpi=600)
