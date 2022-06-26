import numpy

import plotTimeEvol
import plotGSProps
import numpy as np
import plotPolaritonFreqs
import classicalVSQuantum
import ConvergencePlots
import gsProps2Bands

import h5py

def main():
    numpy.set_printoptions(precision=10)

    #calculate bare phonon frequency in case of quadratic coupling
    #omegaPhon = 0.4 * 0.109565595278 / 2 +  2 * np.sqrt(0.2**2 * (0.109565595278 / 2)**2 + 1)
    #print("omegaPhon = {}".format(omegaPhon))

    #plotGSProps.plotGSPropsOnlyPhot()

    filenames = [
        "../data/gsProp2BandsUa50Ub50Uud50Uss50epsA0epsB10gE0.hdf5",
        "../data/gsProp2BandsUa50Ub50Uud50Uss50epsA0epsB10gE1.hdf5",
        "../data/gsProp2BandsUa50Ub50Uud50Uss50epsA0epsB10gE5.hdf5",
        "../data/gsProp2BandsUa50Ub50Uud50Uss50epsA0epsB10gE10.hdf5",
        "../data/gsProp2BandsUa50Ub50Uud50Uss50epsA0epsB10gE20.hdf5"
        #"../data/gsProp2BandsUa50Ub50Uud50Uss50epsA0epsB10gE40.hdf5"
    ]

    gsProps2Bands.dOccAsOfWp(filenames)

#wMinus = plotPolaritonFreqs.calcWMinus(2., 2., 0.5)
    #wPlus = plotPolaritonFreqs.calcWPlus(2., 2., 0.5)

    #exit()

    #print("W+ = {}".format(wPlus))
    #print("W- = {}".format(wMinus))
#
    ##exit()
#
    #plotPolaritonFreqs.plotPolaritonFreqs()

    #plotGSProps.plotGSProps()
    #plotGSProps.plotGSPropsTemp()

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
