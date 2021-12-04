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


def main():
    print("plotting some beautiful time evolution")

    #read in stuff
    file = h5py.File("../data/timeEvolutionResults.hdf5", 'r')
    times = file['times'][()]
    pump = file['pump'][()]
    dOcc = file['dOcc'][()]
    Xpt = file['Xpt'][()]
    XptSqr = file['XptSqr'][()]
    Xph = file['Xph'][()]
    XphSqr = file['XphSqr'][()]


    print("times.shape = {}".format(times.shape))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(times, pump)
    ax.plot(times, Xpt)

    #ax.set_ylim(-1e10, 1e10)

    plt.show()



main()
