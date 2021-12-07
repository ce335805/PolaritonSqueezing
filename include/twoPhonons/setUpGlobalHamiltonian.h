#ifndef PHONONSQUEEZING_SETUPGLOBALHAMILTONIAN_H
#define PHONONSQUEEZING_SETUPGLOBALHAMILTONIAN_H

void setupUncoupledHamiltonian(std::vector<std::complex<double>> &H);

void setupEPhCoupling(std::vector<std::complex<double>> &HEPh);

void setupGlobalH(std::vector<std::complex<double>> &H);

void setupPhPtCoupling(std::vector<std::complex<double>> &HPhPt);

void setupQuadEPhCoupling(std::vector<std::complex<double>> &HEPh);

#endif //PHONONSQUEEZING_SETUPGLOBALHAMILTONIAN_H
