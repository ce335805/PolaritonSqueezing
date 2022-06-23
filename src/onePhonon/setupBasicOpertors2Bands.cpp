#include <complex>
#include <vector>
#include <iostream>

#include "globals.h"
#include "utils.h"
#include "matrixOperations.h"
#include "setupElectronicOperators2Bands.h"
#include "setupBasicOperatorsOnePh.h"

void setupHelectronic2Bands(std::vector<std::complex<double>> &HElectronic) {

  HElectronic = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> HElectronicSmall;
  setupHElectronicSmall2Bands(HElectronicSmall);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          HElectronic[toGlobalMatrixIndexOne(ptInd, ptInd, ph1Ind, ph1Ind, eInd1, eInd2)] = HElectronicSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }

}

void setupDOcc0(std::vector<std::complex<double>> &DOcc) {

  DOcc = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> DOccSmall;
  setupDoubleOccSmall0(DOccSmall);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          DOcc[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)] = DOccSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupDOcc1(std::vector<std::complex<double>> &DOcc) {
  
  DOcc = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> DOccSmall;
  setupDoubleOccSmall1(DOccSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          DOcc[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)] = DOccSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupInterbandHop0(std::vector<std::complex<double>> &hopInter) {
  
  hopInter = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> hopInterSmall;
  setupInterbandHop0Small(hopInterSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          hopInter[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
          = hopInterSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupInterbandHop1(std::vector<std::complex<double>> &hopInter) {
  
  hopInter = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> hopInterSmall;
  setupInterbandHop1Small(hopInterSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          hopInter[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
              = hopInterSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupA0(std::vector<std::complex<double>> &A0) {
  A0 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd1 = 0ul; phInd1 < dimPhonon; ++phInd1) {
      for (ulong phInd2 = 0ul; phInd2 < dimPhonon; ++phInd2) {
        for (ulong eInd = 0ul; eInd < dimElectron; ++eInd) {
          if (phInd2 == phInd1 + 1ul) {
            A0[toGlobalMatrixIndexOne(ptInd, ptInd, phInd1, phInd2, eInd, eInd)] = std::sqrt(double(phInd2));
          }
        }
      }
    }
  }

}

void setupA1(std::vector<std::complex<double>> &A1) {
  A1 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  for (ulong ptInd1 = 0ul; ptInd1 < dimPhoton; ++ptInd1) {
    for (ulong ptInd2 = 0ul; ptInd2 < dimPhoton; ++ptInd2) {
      for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
        for (ulong eInd = 0ul; eInd < dimElectron; ++eInd) {
          if (ptInd2 == ptInd1 + 1ul) {
            A1[toGlobalMatrixIndexOne(ptInd1, ptInd2, phInd, phInd, eInd, eInd)] = std::sqrt(double(ptInd2));
          }
        }
      }
    }
  }

}

