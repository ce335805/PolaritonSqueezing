#include <complex>
#include <vector>
#include <iostream>

#include "globals.h"
#include "utils.h"
#include "matrixOperations.h"
#include "setupElectronicOperatorsSmall.h"
#include "setupBasicOperatorsOnlyPhot.h"

void setupHelectronicOnlyPhot(std::vector<std::complex<double>> &HElectronic) {
  
  HElectronic = std::vector<std::complex<double>>(dimHOnlyPhot * dimHOnlyPhot, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> HElectronicSmall;
  setupHElectronicSmall(HElectronicSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
      for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
        HElectronic[toGlobalMatrixIndexOnlyPhot(ptInd, ptInd, eInd1, eInd2)] = HElectronicSmall[eInd1 * dimElectron +
                                                                                                eInd2];
      }
    }
  }
}

void setupHCouplingOnlyPhot(std::vector<std::complex<double>> &HCoupling) {
  
  HCoupling = std::vector<std::complex<double>>(dimHOnlyPhot * dimHOnlyPhot, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> HCoupingSmall;
  setupHCouplingSmall(HCoupingSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
      for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
        HCoupling[toGlobalMatrixIndexOnlyPhot(ptInd, ptInd, eInd1, eInd2)] = HCoupingSmall[eInd1 * dimElectron +
                                                                                           eInd2];
      }
    }
  }
}

void setupDOccOnlyPhot(std::vector<std::complex<double>> &DOcc) {
  
  DOcc = std::vector<std::complex<double>>(dimHOnlyPhot * dimHOnlyPhot, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> DOccSmall;
  setupDoubleOccSmall(DOccSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
      for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
        DOcc[toGlobalMatrixIndexOnlyPhot(ptInd, ptInd, eInd1, eInd2)] = DOccSmall[eInd1 * dimElectron + eInd2];
      }
    }
  }
}

void setupAOnlyPhot(std::vector<std::complex<double>> &A) {
  A = std::vector<std::complex<double>>(dimHOnlyPhot * dimHOnlyPhot, std::complex<double>(0., 0.));
  
  for (ulong ptInd1 = 0ul; ptInd1 < dimPhoton; ++ptInd1) {
    for (ulong ptInd2 = 0ul; ptInd2 < dimPhoton; ++ptInd2) {
      for (ulong eInd = 0ul; eInd < dimElectron; ++eInd) {
        if (ptInd2 == ptInd1 + 1ul) {
          A[toGlobalMatrixIndexOnlyPhot(ptInd1, ptInd2, eInd, eInd)] = std::sqrt(double(ptInd2));
        }
      }
    }
  }
  
}
