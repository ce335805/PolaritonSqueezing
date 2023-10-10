import numpy

import plotLevels
import plotTimeEvol
import plotGSProps
import numpy as np
import plotPolaritonFreqs
import classicalVSQuantum
import ConvergencePlots
import gsProps2Bands
import timeEvol2Bands
import scipy.constants as consts

import h5py

from interbandPaper import asOfFreq
from interbandPaper import plotTimeEvolCorr
from interbandPaper import plotTimeEvolOcc
from plotPT import plotPTN
from plotPT import plotPTAlpha
from timeEvolConvergence import plotTimeEvolConvergenceNB
from interbandPaper import plotTimeEvolCorrSideBySide

import plotGSProps

def main():
    numpy.set_printoptions(precision=10)

    plotGSProps.plotGSProps()
    plotGSProps.plotGSPropsTemp()

    #asOfFreq()
    #plotTimeEvolCorr()
    #plotTimeEvolCorrSideBySide()
    #plotTimeEvolOcc()

    #plotTimeEvolConvergenceNB()


    #plotPTN()
    #plotPTAlpha()

    #get spectrum lines

    #filename = "../data/spectrum2BandstHop0epsA0epsB100Om100.hdf5"
    #plotLevels.makeLevelPlot(filename)
    #plotLevels.plotEigenstateProperties(filename)



main()
