import numpy

import plotTimeEvol
import plotGSProps
import numpy as np
import plotPolaritonFreqs
import classicalVSQuantum
import ConvergencePlots

import h5py

def main():
    numpy.set_printoptions(precision=1)

    wMinus = plotPolaritonFreqs.calcWMinus(2., 2., 1.6)
    wPlus = plotPolaritonFreqs.calcWPlus(2., 2., 1.6)

    print("W+ = {}".format(wPlus))
    print("W- = {}".format(wMinus))

    #exit()

    #Plot of the polariton frequencies - Fig. 1 b.)
    plotPolaritonFreqs.plotPolaritonFreqs()

    #Plot of the double occupancy in the ground-state - Fig. 1 c.)
    plotGSProps.plotGSProps()

    #Plot of time-evolution for the coupling to the linear density - Fig. 2
    plotTimeEvol.plotLinUSCpSC()

    #Plot of time-evolution for the coupling to the double occupancy - Fig. 3
    plotTimeEvol.plotQuadUSCpSC()


### plots for the appendix

    #Plots for convergence of double occupancy in boson cutoff - Fig. 5
    plotGSProps.plotGSPropsConvergenceLin()
    plotGSProps.plotGSPropsConvergenceQuad()

    #Plots for convergence of time-evolved system in boson cutoff and time-step - Fig. 6
    ConvergencePlots.plotConvInDOccLin()
    ConvergencePlots.plotConvInDOccQuad()

    #Plots for comparison of classical vs cavity driving - Fig. 7
    classicalVSQuantum.plotClassicalVSQuantumLin()
    classicalVSQuantum.plotClassicalVSQuantumQuad()

    #Plots for comparison of 1Phonon vs 2Phonon model - Fig. 8
    plotTimeEvol.plot1Phonvs2Phon()


main()
