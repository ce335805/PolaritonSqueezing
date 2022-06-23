#include <vector>
#include <complex>
#include <iostream>

#include "globals.h"
#include "setupElectronicOperatorsSmall.h"
#include "writeStuffToHdf5.h"


void setupHElectronicSmall2Bands(std::vector<std::complex<double>> &HElectronicSmall) {

  HElectronicSmall = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));

  std::string filename = "./setupH/savedOperators/HU0_1.0_U1_1.0_eps0_0.0_eps1_100.0.hdf5";

  readInComplex2DArray(HElectronicSmall, filename);
  
}


void setupDoubleOccSmall0(std::vector<std::complex<double>> &DOccSmall) {


  DOccSmall = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));

  std::string filename = "./setupH/savedOperators/dOcc0.hdf5";

  readInComplex2DArray(DOccSmall, filename);

}

void setupDoubleOccSmall1(std::vector<std::complex<double>> &DOccSmall) {
  
  
  DOccSmall = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  
  std::string filename = "./setupH/savedOperators/dOcc1.hdf5";
  
  readInComplex2DArray(DOccSmall, filename);
  
}

void setupInterbandHop0Small(std::vector<std::complex<double>> &interbandHop) {
  
  interbandHop = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  
  std::string filename = "./setupH/savedOperators/InterOrb0.hdf5";
  
  readInComplex2DArray(interbandHop, filename);
  
}

void setupInterbandHop1Small(std::vector<std::complex<double>> &interbandHop) {
  
  
  interbandHop = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  
  std::string filename = "./setupH/savedOperators/InterOrb1.hdf5";
  
  readInComplex2DArray(interbandHop, filename);
  
}