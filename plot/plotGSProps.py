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

    fileOnePh = h5py.File("../dataRef/gsPropQuad2PhGPH20NB16.hdf5", 'r')
    wPOnePh = fileOnePh['times'][()]
    dOccOnePh = fileOnePh['dOcc'][()]
    NptOnePh = fileOnePh['Npt'][()]
    NphOnePh = fileOnePh['N1ph'][()]

    fileTwoPh = h5py.File("../dataRef/gsPropQuad2PhGPH50NB16.hdf5", 'r')
    wPTwoPh = fileTwoPh['times'][()]
    dOccTwoPh = fileTwoPh['dOcc'][()]
    NptTwoPh = fileTwoPh['Npt'][()]
    NphTwoPh = fileTwoPh['N1ph'][()]

    fig = plt.figure()
    ax = fig.add_subplot(111)


    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)

    ax.tick_params(direction='inout', length=4, width=.8)

    fig.set_size_inches(2.3, 2.)
    #fig.set_size_inches(5., 5.)

    linewidth = 1.
    markersize = 3
    markeredgewidth = 0.5
    markeredgecolor = 'black'

    ax.plot(wPOnePh / 2. * np.sqrt(2.), dOccOnePh + 0.5, color='olive', label='$X^2 n^2$', marker = 'o', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)
    ax.plot(wPTwoPh / 2. * np.sqrt(2.), dOccTwoPh + 0.5, color='rosybrown', label='$X^2 n$', marker = 's', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)

    ax.plot(wPTwoPh / 2. * np.sqrt(2.), np.ones(len(wPTwoPh)) * 0.109566, color='black', label='$g = 0$', marker = '', linewidth = linewidth)

    #ax.plot(wPTwoPh, NphTwoPh, color='olive', label=r'$\langle X^2 \rangle$ - $X^2 n$')
    #ax.plot(wPTwoPh, NptTwoPh, color='rosybrown', label=r'$\langle X^2 \rangle$ - $X^2 n^2$')

    ax.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize)
    ax.set_xlabel(r"$\rm light{-}matter$ $\rm coupling$ $[\omega_{\rm P} / \omega_{\rm phot}]$", fontsize = fontsize-1)


    ax.set_xlim(0., 3.)

    ax.set_xticks([0, 1., 2., 3.])
    ax.set_xticklabels(["$0$", "$1$", "$2$", "$3$"], fontsize = fontsize)

    ax.set_yticks([0.11, 0.112, 0.114, 0.116])
    ax.set_yticklabels(["$0.11$", "$0.112$", "$0.114$", "$0.116$"], fontsize = fontsize)

    arrow = patches.FancyArrowPatch((.5, 0.1137), (2.5, 0.1127), arrowstyle='->', mutation_scale=10, zorder = 100, linewidth=1.5, color = 'black')
    ax.add_patch(arrow)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=0., linewidth=0., fill=True, pad=0.15)
    #plt.gcf().text(.4, 0.65, "Increased Correlations", fontsize=10, color='black', alpha=1., bbox=boxProps)
    ax.text(.1, 0.114, r"$\rm Increased \, \, \,  Interactions$", fontsize=fontsize - 1, color='black', alpha=1., bbox=boxProps)

    #plt.gcf().text(.9, 0.1, r"$g_1 = g_2 = 0$", fontsize=10, color='black', alpha=1., bbox=boxProps)
    ax.text(1., 0.1099, r"$g = 0$", fontsize=10, color='black', alpha=1., bbox=boxProps)
    ax.text(2.3, 0.1109, r"$X^2 n$", fontsize=10, color='rosybrown', alpha=1., bbox=boxProps)
    ax.text(2., 0.11545, r"$X^2 n_{\uparrow}n_{\downarrow}$", fontsize=10, color='olive', alpha=1., bbox=boxProps)


    #legend = ax.legend(fontsize = fontsize, loc = 'upper left', bbox_to_anchor=(.0, .8), edgecolor = 'black', ncol = 1)
    #legend.get_frame().set_alpha(0.)
    #legend.get_frame().set_boxstyle('Square', pad=0.1)
    #legend.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    plt.show()

    #plt.savefig('dOccAsOfwP.png', format='png', bbox_inches='tight', dpi = 600)


def plotGSPropsConvergence():
    print("plotting GS props")

    fileOnePhN6 = h5py.File("../data/gsPropQuad2PhGPH20NB6.hdf5", 'r')
    wP = fileOnePhN6['times'][()]
    dOccOnePhN6 = fileOnePhN6['dOcc'][()]

    fileOnePhN10 = h5py.File("../data/gsPropQuad2PhGPH20NB10.hdf5", 'r')
    dOccOnePhN10 = fileOnePhN10['dOcc'][()]

    fileOnePhN12 = h5py.File("../data/gsPropQuad2PhGPH20NB12.hdf5", 'r')
    dOccOnePhN12 = fileOnePhN12['dOcc'][()]

    fileOnePhN14 = h5py.File("../data/gsPropQuad2PhGPH20NB14.hdf5", 'r')
    dOccOnePhN14 = fileOnePhN14['dOcc'][()]

    fileOnePhN16 = h5py.File("../data/gsPropQuad2PhGPH20NB16.hdf5", 'r')
    dOccOnePhN16 = fileOnePhN16['dOcc'][()]



    fileTwoPhN6 = h5py.File("../data/gsProp2PhGPH50NB6.hdf5", 'r')
    dOccTwoPhN6 = fileTwoPhN6['dOcc'][()]
    
    fileTwoPhN10 = h5py.File("../data/gsProp2PhGPH50NB10.hdf5", 'r')
    dOccTwoPhN10 = fileTwoPhN10['dOcc'][()]
    
    fileTwoPhN12 = h5py.File("../data/gsProp2PhGPH50NB12.hdf5", 'r')
    dOccTwoPhN12 = fileTwoPhN12['dOcc'][()]
    
    fileTwoPhN14 = h5py.File("../data/gsProp2PhGPH50NB14.hdf5", 'r')
    dOccTwoPhN14 = fileTwoPhN14['dOcc'][()]

    fileTwoPhN16 = h5py.File("../data/gsProp2PhGPH50NB16.hdf5", 'r')
    dOccTwoPhN16 = fileTwoPhN16['dOcc'][()]
    
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex = ax1)


    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(0.5)

    fig.set_size_inches(3.5, 3.)
    #fig.set_size_inches(5., 5.)

    linewidth = 1.5
    markersize = 0
    markeredgewidth = 0.5
    markeredgecolor = 'black'

    import matplotlib

    bone = cm.get_cmap('bone')
    pink = cm.get_cmap('pink')

    ax1.plot(wP / 2., dOccOnePhN6 + 0.5, color=bone(0.8), label=r'$N_{\rm B} = 6$', marker = 'o', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)
    ax1.plot(wP / 2., dOccOnePhN10 + 0.5, color=bone(0.6), label=r'$N_{\rm B} = 10$', marker = 'o', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)
    ax1.plot(wP / 2., dOccOnePhN12 + 0.5, color=bone(0.4), label=r'$N_{\rm B} = 12$', marker = 'o', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)
    ax1.plot(wP / 2., dOccOnePhN14 + 0.5, color=bone(0.2), label=r'$N_{\rm B} = 14$', marker = 'o', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)
    ax1.plot(wP / 2., dOccOnePhN16 + 0.5, color='rosybrown', linestyle = '--', label=r'$N_{\rm B} = 16$', marker = 'o', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)


    ax2.plot(wP / 2., dOccTwoPhN6 + 0.5, color=pink(0.8), label=r'$N_{\rm B} = 6$', marker = 's', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)
    ax2.plot(wP / 2., dOccTwoPhN10 + 0.5, color=pink(0.6), label=r'$N_{\rm B} = 10$', marker = 's', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)
    ax2.plot(wP / 2., dOccTwoPhN12 + 0.5, color=pink(0.4), label=r'$N_{\rm B} = 12$', marker = 's', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)
    ax2.plot(wP / 2., dOccTwoPhN14 + 0.5, color=pink(0.2), label=r'$N_{\rm B} = 14$', marker = 's', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)
    ax2.plot(wP / 2., dOccTwoPhN16 + 0.5, color='olive', linestyle = '--', label=r'$N_{\rm B} = 16$', marker = 's', linewidth = linewidth, markersize = markersize, markeredgewidth = markeredgewidth, markeredgecolor = markeredgecolor)


    #ax.plot(wPTwoPh, NphTwoPh, color='olive', label=r'$\langle X^2 \rangle$ - $X^2 n$')
    #ax.plot(wPTwoPh, NptTwoPh, color='rosybrown', label=r'$\langle X^2 \rangle$ - $X^2 n^2$')

    ax2.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize = fontsize)
    ax2.set_xlabel(r"$\omega_{\rm P} / \omega_{\rm phot}$", fontsize = fontsize)


    ax1.set_xlim(0., 3)

    ax1.set_xticks([])
    ax1.set_xticklabels([])
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax2.set_xticks([0, 1., 2., 3.])
    ax2.set_xticklabels(["$0$", "$1$", "$2$", "$3$"], fontsize = fontsize)

    #ax.set_yticks([0.11, 0.112, 0.114, 0.116])
    #ax.set_yticklabels(["$0.11$", "$0.112$", "$0.114$", "$0.116$"], fontsize = fontsize)

    legend1 = ax1.legend(fontsize = fontsize - 2, loc = 'upper left', bbox_to_anchor=(.0, .75), edgecolor = 'black', ncol = 1)
    legend1.get_frame().set_alpha(0.)
    legend1.get_frame().set_boxstyle('Square', pad=0.1)
    legend1.get_frame().set_linewidth(0.0)

    legend2 = ax2.legend(fontsize = fontsize - 2, loc = 'upper left', bbox_to_anchor=(.0, .75), edgecolor = 'black', ncol = 1)
    legend2.get_frame().set_alpha(0.)
    legend2.get_frame().set_boxstyle('Square', pad=0.1)
    legend2.get_frame().set_linewidth(0.0)


    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    #plt.savefig('gsPropsConvergenceNB.png', format='png', bbox_inches='tight', dpi = 600)


def plotGSPropsConvergenceLin():
    print("plotting GS props")

    filePhN6 = h5py.File("../dataRef/gsPropQuad2PhGPH50NB6.hdf5", 'r')
    wP = filePhN6['times'][()]
    dOccPhN6 = filePhN6['dOcc'][()]
    NphPhN6 = 2. * filePhN6['N1ph'][()]
    NptPhN6 = filePhN6['Npt'][()]

    filePhN10 = h5py.File("../dataRef/gsPropQuad2PhGPH50NB10.hdf5", 'r')
    dOccPhN10 = filePhN10['dOcc'][()]
    NphPhN10 = 2. * filePhN10['N1ph'][()]
    NptPhN10 = filePhN10['Npt'][()]

    filePhN12 = h5py.File("../dataRef/gsPropQuad2PhGPH50NB12.hdf5", 'r')
    dOccPhN12 = filePhN12['dOcc'][()]
    NphPhN12 = 2. * filePhN12['N1ph'][()]
    NptPhN12 = filePhN12['Npt'][()]

    filePhN14 = h5py.File("../dataRef/gsPropQuad2PhGPH50NB14.hdf5", 'r')
    dOccPhN14 = filePhN14['dOcc'][()]
    NphPhN14 = 2. * filePhN14['N1ph'][()]
    NptPhN14 = filePhN14['Npt'][()]

    filePhN16 = h5py.File("../dataRef/gsPropQuad2PhGPH50NB16.hdf5", 'r')
    dOccPhN16 = filePhN16['dOcc'][()]
    NphPhN16 = 2. * filePhN16['N1ph'][()]
    NptPhN16 = filePhN16['Npt'][()]

    fig = plt.figure()
    gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex = ax1)
    ax3 = fig.add_subplot(gs[2], sharex = ax1)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(0.5)
        ax2.spines[axis].set_linewidth(0.5)
        ax3.spines[axis].set_linewidth(0.5)

    fig.set_size_inches(3., 3.)
    # fig.set_size_inches(5., 5.)

    linewidth = 1.5
    markersize = 0
    markeredgewidth = 0.5
    markeredgecolor = 'black'

    import matplotlib

    pink = cm.get_cmap('pink')
    colorFinal = '#184673'

    ax1.plot(wP / 2. * np.sqrt(2.), dOccPhN6 + 0.5, color=pink(0.15), label=r'$N_{\rm B} = 6$', marker='', linewidth=linewidth)
    ax1.plot(wP / 2. * np.sqrt(2.), dOccPhN10 + 0.5, color=pink(0.3), label=r'$N_{\rm B} = 10$', marker='', linewidth=linewidth)
    ax1.plot(wP / 2. * np.sqrt(2.), dOccPhN12 + 0.5, color=pink(0.45), label=r'$N_{\rm B} = 12$', marker='', linewidth=linewidth)
    ax1.plot(wP / 2. * np.sqrt(2.), dOccPhN14 + 0.5, color=pink(0.6), label=r'$N_{\rm B} = 14$', marker='', linewidth=linewidth)
    ax1.plot(wP / 2. * np.sqrt(2.), dOccPhN16 + 0.5, color=colorFinal, linestyle='--', dashes=[4, 4], label=r'$N_{\rm B} = 16$', marker='', linewidth=linewidth)
    ax1.axvline(3., 0., 1., color = 'red', linewidth = 0.5)

    ax2.plot(wP / 2. * np.sqrt(2.), NphPhN6, color=pink(0.15), label=r'$N_{\rm B} = 6$', marker='', linewidth=linewidth)
    ax2.plot(wP / 2. * np.sqrt(2.), NphPhN10, color=pink(0.3), label=r'$N_{\rm B} = 10$', marker='', linewidth=linewidth)
    ax2.plot(wP / 2. * np.sqrt(2.), NphPhN12, color=pink(0.45), label=r'$N_{\rm B} = 12$', marker='', linewidth=linewidth)
    ax2.plot(wP / 2. * np.sqrt(2.), NphPhN14, color=pink(0.6), label=r'$N_{\rm B} = 14$', marker='', linewidth=linewidth)
    ax2.plot(wP / 2. * np.sqrt(2.), NphPhN16, color=colorFinal, linestyle='--', dashes=[4, 4], label=r'$N_{\rm B} = 16$', marker='', linewidth=linewidth)
    ax2.axvline(3., 0., 1., color = 'red', linewidth = 0.5)


    ax3.plot(wP / 2. * np.sqrt(2.), NptPhN6, color=pink(0.15), label=r'$N_{\rm B} = 6$', marker='', linewidth=linewidth)
    ax3.plot(wP / 2. * np.sqrt(2.), NptPhN10, color=pink(0.3), label=r'$N_{\rm B} = 10$', marker='', linewidth=linewidth)
    ax3.plot(wP / 2. * np.sqrt(2.), NptPhN12, color=pink(0.45), label=r'$N_{\rm B} = 12$', marker='', linewidth=linewidth)
    ax3.plot(wP / 2. * np.sqrt(2.), NptPhN14, color=pink(0.6), label=r'$N_{\rm B} = 14$', marker='', linewidth=linewidth)
    ax3.plot(wP / 2. * np.sqrt(2.), NptPhN16, color=colorFinal, linestyle='--', dashes=[4, 4], label=r'$N_{\rm B} = 16$', marker='', linewidth=linewidth)
    ax3.axvline(3., 0., 1., color = 'red', linewidth = 0.5)

    # ax.plot(wPTwoPh, NphTwoPh, color='olive', label=r'$\langle X^2 \rangle$ - $X^2 n$')
    # ax.plot(wPTwoPh, NptTwoPh, color='rosybrown', label=r'$\langle X^2 \rangle$ - $X^2 n^2$')

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize=fontsize, labelpad=22.)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize=fontsize, labelpad=22.)


    ax3.set_xlabel(r"$\omega_{\rm P} / \omega_{\rm phot}$", fontsize=fontsize)

    ax1.set_xlim(0., 6.)
    ax2.set_xlim(0., 6.)
    ax3.set_xlim(0., 6.)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax3.set_xticks([0, 2., 4., 6.])
    ax3.set_xticklabels(["$0$", "$2$", "$4$", "$6$"], fontsize=fontsize)

    ax1.set_yticks([0.110, 0.111])
    ax1.set_yticklabels(["$0.110$", "$0.111$"], fontsize = fontsize)

    ax2.set_ylim(-0.2, 2.1)
    ax2.set_yticks([0., 2.])
    ax2.set_yticklabels(["$0$", "$2$"], fontsize = fontsize)

    ax3.set_ylim(-0.2, 1.3)
    ax3.set_yticks([0., 1.])
    ax3.set_yticklabels(["$0$", "$1$"], fontsize = fontsize)

    arrow = patches.FancyArrowPatch((2.99, 0.1098), (2.2, 0.1098), arrowstyle='->', mutation_scale=10, zorder = 100, linewidth=1., color = 'black')
    ax1.add_patch(arrow)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax1.text(0.4, 0.11, r"$\rm Values$ $\rm in$ $\rm paper$", fontsize=8, color='black', alpha=1., bbox=boxProps)

    legend2 = ax1.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(.0, 1.6), edgecolor='black', ncol=2)
    legend2.get_frame().set_alpha(0.)
    legend2.get_frame().set_boxstyle('Square', pad=0.1)
    legend2.get_frame().set_linewidth(0.0)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax1.text(-0.5, 0.1115, r"$\rm{\textbf{a.)}}$", fontsize=10, color='black', alpha=1., bbox=boxProps)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    #plt.savefig('gsPropsConvergenceNBLin.png', format='png', bbox_inches='tight', dpi = 600)


def plotGSPropsConvergenceQuad():
    print("plotting GS props")

    filePhN6 = h5py.File("../dataRef/gsPropQuad2PhGPH20NB6.hdf5", 'r')
    wP = filePhN6['times'][()]
    dOccPhN6 = filePhN6['dOcc'][()]
    NphPhN6 = 2. * filePhN6['N1ph'][()]
    NptPhN6 = filePhN6['Npt'][()]

    filePhN10 = h5py.File("../dataRef/gsPropQuad2PhGPH20NB10.hdf5", 'r')
    dOccPhN10 = filePhN10['dOcc'][()]
    NphPhN10 = 2. * filePhN10['N1ph'][()]
    NptPhN10 = filePhN10['Npt'][()]

    filePhN12 = h5py.File("../dataRef/gsPropQuad2PhGPH20NB12.hdf5", 'r')
    dOccPhN12 = filePhN12['dOcc'][()]
    NphPhN12 = 2. * filePhN12['N1ph'][()]
    NptPhN12 = filePhN12['Npt'][()]

    filePhN14 = h5py.File("../dataRef/gsPropQuad2PhGPH20NB14.hdf5", 'r')
    dOccPhN14 = filePhN14['dOcc'][()]
    NphPhN14 = 2. * filePhN14['N1ph'][()]
    NptPhN14 = filePhN14['Npt'][()]

    filePhN16 = h5py.File("../dataRef/gsPropQuad2PhGPH20NB16.hdf5", 'r')
    dOccPhN16 = filePhN16['dOcc'][()]
    NphPhN16 = 2. * filePhN16['N1ph'][()]
    NptPhN16 = filePhN16['Npt'][()]

    fig = plt.figure()
    gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex = ax1)
    ax3 = fig.add_subplot(gs[2], sharex = ax1)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(0.5)
        ax2.spines[axis].set_linewidth(0.5)
        ax3.spines[axis].set_linewidth(0.5)

    fig.set_size_inches(3., 2.62)
    # fig.set_size_inches(5., 5.)

    linewidth = 1.5
    markersize = 0
    markeredgewidth = 0.5
    markeredgecolor = 'black'

    import matplotlib

    pink = cm.get_cmap('pink')
    colorFinal = '#184673'

    ax1.plot(wP / 2. * np.sqrt(2.), dOccPhN6 + 0.5, color=pink(0.15), label=r'$N_{\rm B} = 6$', marker='', linewidth=linewidth)
    ax1.plot(wP / 2. * np.sqrt(2.), dOccPhN10 + 0.5, color=pink(0.3), label=r'$N_{\rm B} = 10$', marker='', linewidth=linewidth)
    ax1.plot(wP / 2. * np.sqrt(2.), dOccPhN12 + 0.5, color=pink(0.45), label=r'$N_{\rm B} = 12$', marker='', linewidth=linewidth)
    ax1.plot(wP / 2. * np.sqrt(2.), dOccPhN14 + 0.5, color=pink(0.6), label=r'$N_{\rm B} = 14$', marker='', linewidth=linewidth)
    ax1.plot(wP / 2. * np.sqrt(2.), dOccPhN16 + 0.5, color=colorFinal, linestyle='--', dashes=[4, 4], label=r'$N_{\rm B} = 16$', marker='', linewidth=linewidth)
    ax1.axvline(3., 0., 1., color = 'red', linewidth = 0.5)

    ax2.plot(wP / 2. * np.sqrt(2.), NphPhN6, color=pink(0.15), label=r'$N_{\rm B} = 6$', marker='', linewidth=linewidth)
    ax2.plot(wP / 2. * np.sqrt(2.), NphPhN10, color=pink(0.3), label=r'$N_{\rm B} = 10$', marker='', linewidth=linewidth)
    ax2.plot(wP / 2. * np.sqrt(2.), NphPhN12, color=pink(0.45), label=r'$N_{\rm B} = 12$', marker='', linewidth=linewidth)
    ax2.plot(wP / 2. * np.sqrt(2.), NphPhN14, color=pink(0.6), label=r'$N_{\rm B} = 14$', marker='', linewidth=linewidth)
    ax2.plot(wP / 2. * np.sqrt(2.), NphPhN16, color=colorFinal, linestyle='--', dashes=[4, 4], label=r'$N_{\rm B} = 16$', marker='', linewidth=linewidth)
    ax2.axvline(3., 0., 1., color = 'red', linewidth = 0.5)


    ax3.plot(wP / 2. * np.sqrt(2.), NptPhN6, color=pink(0.15), label=r'$N_{\rm B} = 6$', marker='', linewidth=linewidth)
    ax3.plot(wP / 2. * np.sqrt(2.), NptPhN10, color=pink(0.3), label=r'$N_{\rm B} = 10$', marker='', linewidth=linewidth)
    ax3.plot(wP / 2. * np.sqrt(2.), NptPhN12, color=pink(0.45), label=r'$N_{\rm B} = 12$', marker='', linewidth=linewidth)
    ax3.plot(wP / 2. * np.sqrt(2.), NptPhN14, color=pink(0.6), label=r'$N_{\rm B} = 14$', marker='', linewidth=linewidth)
    ax3.plot(wP / 2. * np.sqrt(2.), NptPhN16, color=colorFinal, linestyle='--', dashes=[4, 4], label=r'$N_{\rm B} = 16$', marker='', linewidth=linewidth)
    ax3.axvline(3., 0., 1., color = 'red', linewidth = 0.5)

    ax1.set_ylabel(r"$\sum_i \langle n_{i, \uparrow} n_{i, \downarrow} \rangle$", fontsize=fontsize)
    ax2.set_ylabel(r"$N_{\rm phon}$", fontsize=fontsize, labelpad=22.)
    ax3.set_ylabel(r"$N_{\rm phot}$", fontsize=fontsize, labelpad=22.)


    ax3.set_xlabel(r"$\omega_{\rm P} / \omega_{\rm phot}$", fontsize=fontsize)

    ax1.set_xlim(0., 6.)
    ax2.set_xlim(0., 6.)
    ax3.set_xlim(0., 6.)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax3.set_xticks([0, 2., 4., 6.])
    ax3.set_xticklabels(["$0$", "$2$", "$4$", "$6$"], fontsize=fontsize)

    ax1.set_yticks([0.114, 0.115, 0.116])
    ax1.set_yticklabels(["$0.114$", "$0.115$", "$0.116$"], fontsize = fontsize)

    ax2.set_ylim(-0.2, 2.1)
    ax2.set_yticks([0., 2.])
    ax2.set_yticklabels(["$0$", "$2$"], fontsize = fontsize)

    ax3.set_ylim(-0.2, 1.3)
    ax3.set_yticks([0., 1.])
    ax3.set_yticklabels(["$0$", "$1$"], fontsize = fontsize)


    arrow = patches.FancyArrowPatch((2.99, 0.114), (2.2, 0.114), arrowstyle='->', mutation_scale=10, zorder = 100, linewidth=1., color = 'black')
    ax1.add_patch(arrow)

    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax1.text(0.4, 0.1143, r"$\rm Values$ $\rm in$ $\rm paper$", fontsize=8, color='black', alpha=1., bbox=boxProps)

    #legend1 = ax1.legend(fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(.0, .7), edgecolor='black', ncol=1)
    #legend1.get_frame().set_alpha(0.)
    #legend1.get_frame().set_boxstyle('Square', pad=0.1)
    #legend1.get_frame().set_linewidth(0.0)


    boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
    ax1.text(-0., 0.1167, r"$\rm{\textbf{b.)}}$", fontsize=10, color='black', alpha=1., bbox=boxProps)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.0)
    plt.show()

    #plt.savefig('gsPropsConvergenceNBQuad.png', format='png', bbox_inches='tight', dpi = 600)
