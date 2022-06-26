#include <complex>
#include <vector>
#include <iostream>

#include "globals.h"
#include "utils.h"
#include "matrixOperations.h"
#include "setupElectronicOperatorsSmall.h"
#include "setupBasicOperatorsOnePh.h"

void setupHelectronicOnePh(std::vector<std::complex<double>> &HElectronic) {

  HElectronic = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> HElectronicSmall;
  setupHElectronicSmall(HElectronicSmall);

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

void setupDOccOnePh(std::vector<std::complex<double>> &DOcc) {

  DOcc = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> DOccSmall;
  setupDoubleOccSmall(DOccSmall);

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

void setupDOccOnePhNoPHS(std::vector<std::complex<double>> &DOcc) {

  DOcc = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> DOccSmall;
  setupDoubleOccSmallNoPHS(DOccSmall);

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

void setupOccSiteIOnePh(std::vector<std::complex<double>> &OccSiteI, const ulong site) {

  OccSiteI = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> occSmall;
  setupOccSiteISmall(occSmall, site);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          OccSiteI[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)] = occSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }

}

void setupB(std::vector<std::complex<double>> &B) {
  B = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd1 = 0ul; phInd1 < dimPhonon; ++phInd1) {
      for (ulong phInd2 = 0ul; phInd2 < dimPhonon; ++phInd2) {
        for (ulong eInd = 0ul; eInd < dimElectron; ++eInd) {
          if (phInd2 == phInd1 + 1ul) {
            B[toGlobalMatrixIndexOne(ptInd, ptInd, phInd1, phInd2, eInd, eInd)] = std::sqrt(double(phInd2));
          }
        }
      }
    }
  }

}

void setupAOnePh(std::vector<std::complex<double>> &A) {
  A = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  for (ulong ptInd1 = 0ul; ptInd1 < dimPhoton; ++ptInd1) {
    for (ulong ptInd2 = 0ul; ptInd2 < dimPhoton; ++ptInd2) {
      for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
        for (ulong eInd = 0ul; eInd < dimElectron; ++eInd) {
          if (ptInd2 == ptInd1 + 1ul) {
            A[toGlobalMatrixIndexOne(ptInd1, ptInd2, phInd, phInd, eInd, eInd)] = std::sqrt(double(ptInd2));
          }
        }
      }
    }
  }

}

