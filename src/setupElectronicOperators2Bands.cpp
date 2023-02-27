#include <vector>
#include <complex>
#include <iostream>

#include "globals.h"
#include "setupElectronicOperatorsSmall.h"
#include "writeStuffToHdf5.h"

#include "matrixOperations.h"

void setupHElectronicSmall2Bands(std::vector<std::complex<double>> &HElectronicSmall) {
  
  HElectronicSmall = std::vector<std::complex<double>> (dimElectron * dimElectron,std::complex<double> (0., 0.));
  
  std::string filename;
  
  //std::string filename = "./setupH/savedOperators/HU0_3.0_U1_3.0_eps0_0.0_eps1_1.0.hdf5";
  //readInComplex2DArray(HElectronicSmall, filename);
  
  std::vector<std::complex<double>> HKin (dimElectron * dimElectron, std::complex<double> (0., 0.));
  filename = "./setupH/savedOperators/HKin.hdf5";
  readInComplex2DArray(HKin, filename);
  
  std::vector<std::complex<double>> HIntBand0 (dimElectron * dimElectron, std::complex<double> (0., 0.));
  filename = "./setupH/savedOperators/dOcc0.hdf5";
  readInComplex2DArray(HIntBand0, filename);
  std::vector<std::complex<double>> HIntBand1 (dimElectron * dimElectron, std::complex<double> (0., 0.));
  filename = "./setupH/savedOperators/dOcc1.hdf5";
  readInComplex2DArray(HIntBand1, filename);
  
  std::vector<std::complex<double>> nSite0 (dimElectron * dimElectron, std::complex<double> (0., 0.));
  filename = "./setupH/savedOperators/onsitePot0.hdf5";
  readInComplex2DArray(nSite0, filename);
  std::vector<std::complex<double>> nSite1 (dimElectron * dimElectron, std::complex<double> (0., 0.));
  filename = "./setupH/savedOperators/onsitePot1.hdf5";
  readInComplex2DArray(nSite1, filename);
  
  std::vector<std::complex<double>> intOrbUpDn (dimElectron * dimElectron, std::complex<double> (0., 0.));
  filename = "./setupH/savedOperators/intOrbUpDn.hdf5";
  readInComplex2DArray(intOrbUpDn, filename);
  
  std::vector<std::complex<double>> intOrbSigSig (dimElectron * dimElectron, std::complex<double> (0., 0.));
  filename = "./setupH/savedOperators/intOrbSigSig.hdf5";
  readInComplex2DArray(intOrbSigSig, filename);
  
  std::vector<std::complex<double>> hPot (dimElectron * dimElectron, std::complex<double> (0., 0.));
  std::vector<std::complex<double>> hInt (dimElectron * dimElectron, std::complex<double> (0., 0.));
  
  addMatricies(nSite0, eps0, nSite1, eps1, hPot);
  addMatricies(HIntBand0, U, HIntBand1, U1, hInt);
  addMatricies(hInt, 1., intOrbUpDn, uUpDn, hInt);
  addMatricies(hInt, 1., intOrbSigSig, uSigSig, hInt);
  
  addMatricies(hInt, hPot, HElectronicSmall);
  addMatricies(HElectronicSmall, 1., HKin, -tHop, HElectronicSmall);
  
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

void setupInterOrbUpDnSmall(std::vector<std::complex<double>> &interOrbUpDn) {
  
  interOrbUpDn = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  
  std::string filename = "./setupH/savedOperators/intOrbUpDn.hdf5";
  
  readInComplex2DArray(interOrbUpDn, filename);
}

void setupInterOrbSigSigSmall(std::vector<std::complex<double>> &interOrbSigSig) {
  
  interOrbSigSig = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  
  std::string filename = "./setupH/savedOperators/intOrbSigSig.hdf5";
  
  readInComplex2DArray(interOrbSigSig, filename);
}

void setupN0Small(std::vector<std::complex<double>> &N0) {
  
  N0 = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  
  std::string filename = "./setupH/savedOperators/n0.hdf5";
  
  readInComplex2DArray(N0, filename);
}

void setupN1Small(std::vector<std::complex<double>> &N1) {
  
  N1 = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  
  std::string filename = "./setupH/savedOperators/n1.hdf5";
  
  readInComplex2DArray(N1, filename);
}

void setupNC0Small(std::vector<std::complex<double>> &NC0) {
  NC0 = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  std::string filename = "./setupH/savedOperators/nc0.hdf5";
  readInComplex2DArray(NC0, filename);
}

void setupND0Small(std::vector<std::complex<double>> &ND0) {
  ND0 = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  std::string filename = "./setupH/savedOperators/nd0.hdf5";
  readInComplex2DArray(ND0, filename);
}

void setupNC1Small(std::vector<std::complex<double>> &NC1) {
  NC1 = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  std::string filename = "./setupH/savedOperators/nc1.hdf5";
  readInComplex2DArray(NC1, filename);
}

void setupND1Small(std::vector<std::complex<double>> &ND1) {
  ND1 = std::vector<std::complex<double>> (dimElectron * dimElectron, std::complex<double> (0., 0.));
  std::string filename = "./setupH/savedOperators/nd1.hdf5";
  readInComplex2DArray(ND1, filename);
}