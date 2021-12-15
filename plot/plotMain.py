import numpy

import plotTimeEvol
import plotGSProps
import numpy as np
import plotPolaritonFreqs

import h5py

def main():

    #plotPolaritonFreqs.plotPolaritonFreqs()

    numpy.set_printoptions(precision=1)

    #print(np.diag(H1)[190:200])
    #print("\n \n \n")
    #print(np.diag(H2)[190:200])

    #plotTimeEvol.plotQCCompareNSqr()
    #plotTimeEvol.plotAsOfWP()
    #plotTimeEvol.plotAsOfWPLin()
    plotGSProps.plotGSProps()


main()
