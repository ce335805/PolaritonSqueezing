#include <complex>
#include <vector>
#include <iostream>

#include "globals.h"
#include "utils.h"
#include "matrixOperations.h"
#include "include/twoPhonons/setupBasicOperators.h"
#include "setupElectronicOperatorsSmall.h"


void setupHelectronic(std::vector<std::complex<double>> &HElectronic) {

  HElectronic = std::vector<std::complex<double>>(dimHTwoPh * dimHTwoPh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> HElectronicSmall;
  setupHElectronicSmall(HElectronicSmall);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
        for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
          for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
            HElectronic[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd1, eInd2)] = HElectronicSmall[eInd1 * dimElectron + eInd2];
          }
        }
      }
    }
  }

}



void setupDOcc(std::vector<std::complex<double>> &DOcc) {

  DOcc = std::vector<std::complex<double>>(dimHTwoPh * dimHTwoPh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> DOccSmall;
  setupDoubleOccSmall(DOccSmall);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
        for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
          for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
            DOcc[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd1, eInd2)] = DOccSmall[eInd1 * dimElectron + eInd2];
          }
        }
      }
    }
  }

}

void setupDOccNoPHS(std::vector<std::complex<double>> &DOcc) {

  DOcc = std::vector<std::complex<double>>(dimHTwoPh * dimHTwoPh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> DOccSmall;
  setupDoubleOccSmallNoPHS(DOccSmall);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
        for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
          for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
            DOcc[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd1, eInd2)] = DOccSmall[eInd1 * dimElectron + eInd2];
          }
        }
      }
    }
  }

}


void setupTotalSpin(std::vector<std::complex<double>> &totalSpin) {

  totalSpin = std::vector<std::complex<double>>(dimHTwoPh * dimHTwoPh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> totalSpinSmall;
  setupTotalSpinSmall(totalSpinSmall);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
        for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
          for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
            totalSpin[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd1, eInd2)] = totalSpinSmall[eInd1 * dimElectron + eInd2];
          }
        }
      }
    }
  }

}

void setupDOccSiteI(std::vector<std::complex<double>> &DOcc, const ulong site) {

  DOcc = std::vector<std::complex<double>>(dimHTwoPh * dimHTwoPh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> DOccSmall;
  setupDOccSiteISmall(DOccSmall, site);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
        for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
          for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
            DOcc[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd1, eInd2)] = DOccSmall[eInd1 * dimElectron + eInd2];
          }
        }
      }
    }
  }

}


void setupOccSiteI(std::vector<std::complex<double>> &OccSiteI, const ulong site) {

  OccSiteI = std::vector<std::complex<double>>(dimHTwoPh * dimHTwoPh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> occSmall;
  setupOccSiteISmall(occSmall, site);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
        for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
          for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
            OccSiteI[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd1, eInd2)] = occSmall[eInd1 * dimElectron + eInd2];
          }
        }
      }
    }
  }

}

void setupB1(std::vector<std::complex<double>> &B1) {
  B1 = std::vector<std::complex<double>>(dimHTwoPh * dimHTwoPh, std::complex<double>(0., 0.));

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind1 = 0ul; ph1Ind1 < dimPhonon; ++ph1Ind1) {
        for (ulong ph1Ind2 = 0ul; ph1Ind2 < dimPhonon; ++ph1Ind2) {
          for (ulong eInd = 0ul; eInd < dimElectron; ++eInd) {
            if (ph1Ind2 == ph1Ind1 + 1ul) {
              B1[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind, ph2Ind, ph1Ind1, ph1Ind2, eInd, eInd)] = std::sqrt(double(ph1Ind2));
            }
          }
        }
      }
    }
  }

}


void setupB2(std::vector<std::complex<double>> &B2) {
  B2 = std::vector<std::complex<double>>(dimHTwoPh * dimHTwoPh, std::complex<double>(0., 0.));

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind1 = 0ul; ph2Ind1 < dimPhonon; ++ph2Ind1) {
      for (ulong ph2Ind2 = 0ul; ph2Ind2 < dimPhonon; ++ph2Ind2) {
        for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
          for (ulong eInd = 0ul; eInd < dimElectron; ++eInd) {
            if (ph2Ind2 == ph2Ind1 + 1ul) {
              B2[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind1, ph2Ind2, ph1Ind, ph1Ind, eInd, eInd)] = std::sqrt(double(ph2Ind2));
            }
          }
        }
      }
    }
  }

}

void setupA(std::vector<std::complex<double>> &A) {
  A = std::vector<std::complex<double>>(dimHTwoPh * dimHTwoPh, std::complex<double>(0., 0.));

  for (ulong ptInd1 = 0ul; ptInd1 < dimPhoton; ++ptInd1) {
    for (ulong ptInd2 = 0ul; ptInd2 < dimPhoton; ++ptInd2) {
      for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
        for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
          for (ulong eInd = 0ul; eInd < dimElectron; ++eInd) {
            if (ptInd2 == ptInd1 + 1ul) {
              A[toGlobalMatrixIndex(ptInd1, ptInd2, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd, eInd)] = std::sqrt(double(ptInd2));
            }
          }
        }
      }
    }
  }

}

