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
}

void setupDoubleOccSmall(std::vector<std::complex<double>> &DOccSmall) {

  DOccSmall = std::vector<std::complex<double>>(
          {0.5, 0., 0., 0.,
           0., -0.5, 0., 0.,
           0., 0., -0.5, 0.,
           0, 0., 0., 0.5}
  );
}

void setupDoubleOccSmallNoPHS(std::vector<std::complex<double>> &DOccSmall) {

  DOccSmall = std::vector<std::complex<double>>(
          {1., 0., 0., 0.,
           0., 0., 0., 0.,
           0., 0., 0., 0.,
           0, 0., 0., 1.}
  );
}

void setupDOccSiteISmall(std::vector<std::complex<double>> &dOccSmallSiteI, const ulong site) {


  if (site == 0ul) {
    dOccSmallSiteI = std::vector<std::complex<double>>(
            {1., 0., 0., 0.,
             0., 0., 0., 0.,
             0., 0., 0., 0.,
             0, 0., 0., 0.}
    );
  } else if (site == 1ul) {
    dOccSmallSiteI = std::vector<std::complex<double>>(
            {0., 0., 0., 0.,
             0., 0., 0., 0.,
             0., 0., 0., 0.,
             0, 0., 0., 1.}
    );
  } else {
    std::cout << "Hubbard dimer only has 2 sites!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << '\n';
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
}

