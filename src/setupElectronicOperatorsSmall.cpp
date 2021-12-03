#include <vector>
#include <complex>
#include <iostream>

#include "globals.h"
#include "setupElectronicOperatorsSmall.h"
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

