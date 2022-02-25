import numpy

import plotTimeEvol
import plotGSProps
import numpy as np
import plotPolaritonFreqs
import classicalVSQuantum
import ConvergencePlots

import h5py

def main():
    numpy.set_printoptions(precision=10)

    #wMinus = plotPolaritonFreqs.calcWMinus(2., 2., 1.6)
    #wPlus = plotPolaritonFreqs.calcWPlus(2., 2., 1.6)
#
    #print("W+ = {}".format(wPlus))
    #print("W- = {}".format(wMinus))
#
    ##exit()
#
    #plotPolaritonFreqs.plotPolaritonFreqs()

    #plotGSProps.plotGSProps()
    plotGSProps.plotGSPropsTemp()

    #plotTimeEvol.plotQuadUSCpSC()
    #plotTimeEvol.plotLinUSCpSC()

    #ConvergencePlots.plotConvInDOccLin()
    #ConvergencePlots.plotConvInDOccQuad()

    #plotTimeEvol.plot1Phonvs2Phon()


    #plotGSProps.plotGSPropsConvergenceLin()
    #plotGSProps.plotGSPropsConvergenceQuad()

    #classicalVSQuantum.plotClassicalVSQuantumQuad()
    #classicalVSQuantum.plotClassicalVSQuantumLin()

main()
