#include <complex>
#include <vector>
#include <iostream>

#include "globals.h"
#include "utils.h"
#include "matrixOperations.h"
#include "include/twoPhonons/setupBasicOperators.h"


void setupHelectronic(std::vector<std::complex<double>> &HElectronic) {

  HElectronic = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> HElectronicSmall;
  setupHElectronicSmall(HElectronicSmall);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
        for (ulong eInd1 = 0ul; eInd1 < 4ul; ++eInd1) {
          for (ulong eInd2 = 0ul; eInd2 < 4ul; ++eInd2) {
            HElectronic[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd1, eInd2)] = HElectronicSmall[eInd1 * 4ul + eInd2];
          }
        }
      }
    }
  }

}

void setupHElectronicSmall(std::vector<std::complex<double>> &HElectronicSmall) {

  HElectronicSmall = std::vector<std::complex<double>>(
          {U / 2., -tHop, tHop, 0.,
           -tHop, -U / 2., 0., -tHop,
           tHop, 0., -U / 2., tHop,
           0, -tHop, tHop, U / 2.}
  );

  //std::cout << "electronic Hamiltonian: " << '\n';
  //for(ulong ind1 = 0ul; ind1 < 4ul; ++ind1){
  //  for(ulong ind2 = 0ul; ind2 < 4ul; ++ind2) {
  //    std::cout << HElectronicSmall[ind1 * 4ul + ind2] << " ";
  //  }
  //  std::cout << '\n';
  //}
}

void setupDoubleOccSmall(std::vector<std::complex<double>> &DOccSmall) {

  DOccSmall = std::vector<std::complex<double>>(
          {0.5, 0., 0., 0.,
           0., -0.5, 0., 0.,
           0., 0., -0.5, 0.,
           0, 0., 0., 0.5}
  );
}

void setupDOcc(std::vector<std::complex<double>> &DOcc) {

  DOcc = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> DOccSmall;
  setupDoubleOccSmall(DOccSmall);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
        for (ulong eInd1 = 0ul; eInd1 < 4ul; ++eInd1) {
          for (ulong eInd2 = 0ul; eInd2 < 4ul; ++eInd2) {
            DOcc[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd1, eInd2)] = DOccSmall[eInd1 * 4ul + eInd2];
          }
        }
      }
    }
  }

}

void setupOccSiteISmall(std::vector<std::complex<double>> &occSmallSiteI, const ulong site) {


  if (site == 0ul) {
    occSmallSiteI = std::vector<std::complex<double>>(
            {2., 0., 0., 0.,
             0., 1., 0., 0.,
             0., 0., 1., 0.,
             0, 0., 0., 0.}
    );
  } else if (site == 1ul) {
    occSmallSiteI = std::vector<std::complex<double>>(
            {0., 0., 0., 0.,
             0., 1., 0., 0.,
             0., 0., 1., 0.,
             0, 0., 0., 2.}
    );
  } else {
    std::cout << "Hubbard dimer only has 2 sites!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << '\n';
  }

  //for (ulong eInd1 = 0ul; eInd1 < 4ul; ++eInd1) {
  //  for (ulong eInd2 = 0ul; eInd2 < 4ul; ++eInd2) {
  //    if (eInd1 == eInd2) {
  //      if (site == 1ul || site == 2ul) {
  //        occSmallSiteI[eInd1 * 4ul + eInd2] = 1.;
  //      }
  //      if (site == 0ul) {
  //        if (eInd1 == 0ul) {
  //          occSmallSiteI[eInd1 * 4ul + eInd2] = 2.;
  //        }
  //      } else if (site == 1ul) {
  //        if (eInd1 == 3ul)
  //          occSmallSiteI[eInd1 * 4ul + eInd2] = 2.;
  //      }
  //    }
  //  }
  //}
}


void setupOccSiteI(std::vector<std::complex<double>> &OccSiteI, const ulong site) {

  OccSiteI = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> occSmall;
  setupOccSiteISmall(occSmall, site);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
        for (ulong eInd1 = 0ul; eInd1 < 4ul; ++eInd1) {
          for (ulong eInd2 = 0ul; eInd2 < 4ul; ++eInd2) {
            OccSiteI[toGlobalMatrixIndex(ptInd, ptInd, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd1, eInd2)] = occSmall[eInd1 * 4ul + eInd2];
          }
        }
      }
    }
  }

}

void setupB1(std::vector<std::complex<double>> &B1) {
  B1 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
      for (ulong ph1Ind1 = 0ul; ph1Ind1 < dimPhonon; ++ph1Ind1) {
        for (ulong ph1Ind2 = 0ul; ph1Ind2 < dimPhonon; ++ph1Ind2) {
          for (ulong eInd = 0ul; eInd < 4ul; ++eInd) {
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
  B2 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph2Ind1 = 0ul; ph2Ind1 < dimPhonon; ++ph2Ind1) {
      for (ulong ph2Ind2 = 0ul; ph2Ind2 < dimPhonon; ++ph2Ind2) {
        for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
          for (ulong eInd = 0ul; eInd < 4ul; ++eInd) {
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
  A = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  for (ulong ptInd1 = 0ul; ptInd1 < dimPhoton; ++ptInd1) {
    for (ulong ptInd2 = 0ul; ptInd2 < dimPhoton; ++ptInd2) {
      for (ulong ph2Ind = 0ul; ph2Ind < dimPhonon; ++ph2Ind) {
        for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
          for (ulong eInd = 0ul; eInd < 4ul; ++eInd) {
            if (ptInd2 == ptInd1 + 1ul) {
              A[toGlobalMatrixIndex(ptInd1, ptInd2, ph2Ind, ph2Ind, ph1Ind, ph1Ind, eInd, eInd)] = std::sqrt(double(ptInd2));
            }
          }
        }
      }
    }
  }

}

