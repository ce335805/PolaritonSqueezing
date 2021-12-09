import numpy

import plotTimeEvol
import plotGSProps
import numpy as np

import h5py

def main():

    numpy.set_printoptions(precision=1)

    file = h5py.File("../data/HTest1.hdf5", 'r')
    H1 = file['HTest'][()]
    file = h5py.File("../data/HTest2.hdf5", 'r')
    H2 = file['HTest'][()]

    #print(H1)
    #print(H2)

    HDiff = H1 - H2

    HPrint = H1

    #for ind1 in range(len(H1)):
    #    for ind2 in range(len(H1)):
    #        if(HPrint[ind1, ind2].real < 0):
    #            print("{0:0.3f} ".format(HPrint[ind1, ind2]), end='')
    #        else:
    #            print(" {0:0.3f} ".format(HPrint[ind1, ind2]), end='')
    #    print('')

    for ind1 in range(len(H1)):
        for ind2 in range(len(H1)):
            if(abs(H1[ind1, ind2] - H2[ind1, ind2]) > 1e-12):
                print("[{}, {}]".format(ind1, ind2))




    #print(np.diag(H1)[190:200])
    #print("\n \n \n")
    #print(np.diag(H2)[190:200])

    plotTimeEvol.plotTimeEvol()
    #plotGSProps.plotGSProps()


main()
