import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
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


def plotGSProps():
    print("plotting GS props")

    fileOnePh = h5py.File("../data/gsPropQuad2PhGPH20NB6.hdf5", 'r')
    wPOnePh = fileOnePh['times'][()]
    dOccOnePh = fileOnePh['dOcc'][()]
    NptOnePh = fileOnePh['Npt'][()]
    NphOnePh = fileOnePh['N1ph'][()]

    fileTwoPh = h5py.File("../data/gsProp2PhGPH50NB6.hdf5", 'r')
    wPTwoPh = fileTwoPh['times'][()]
    dOccTwoPh = fileTwoPh['dOcc'][()]
    NptTwoPh = fileTwoPh['Npt'][()]
    NphTwoPh = fileTwoPh['N1ph'][()]

    fig = plt.figure()
    ax = fig.add_subplot(111)


    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)

    fig.set_size_inches(2.5, 2.)
    #fig.set_size_inches(5., 5.)

    linewidth = 1.
    markersize = 3
    markeredgewidth = 0.5
    markeredgecolor = 'black'

    ax.plot(wPOnePh / 2., dOccOnePh + 0.5, color='olive', label='$X^2 n^2$', marker = 'o', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)
    ax.plot(wPTwoPh / 2., dOccTwoPh + 0.5, color='rosybrown', label='$X^2 n$', marker = 's', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)

    ax.plot(wPTwoPh / 2., np.ones(len(wPTwoPh)) * 0.109566, color='black', label='$g = 0$', marker = '', linewidth = linewidth)

    #ax.plot(wPTwoPh, NphTwoPh, color='olive', label=r'$\langle X^2 \rangle$ - $X^2 n$')
    #ax.plot(wPTwoPh, NptTwoPh, color='rosybrown', label=r'$\langle X^2 \rangle$ - $X^2 n^2$')

    ax.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize)
    ax.set_xlabel(r"$\omega_{\rm P} / \omega_{\rm phot}$", fontsize = fontsize)


    ax.set_xlim(0., 3)

    ax.set_xticks([0, 1., 2., 3.])
    ax.set_xticklabels(["$0$", "$1$", "$2$", "$3$"], fontsize = fontsize)

    ax.set_yticks([0.11, 0.112, 0.114, 0.116])
    ax.set_yticklabels(["$0.11$", "$0.112$", "$0.114$", "$0.116$"], fontsize = fontsize)

    legend = ax.legend(fontsize = fontsize, loc = 'upper left', bbox_to_anchor=(.0, .8), edgecolor = 'black', ncol = 1)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    #plt.show()

    plt.savefig('dOccAsOfwP.png', format='png', bbox_inches='tight', dpi = 600)
