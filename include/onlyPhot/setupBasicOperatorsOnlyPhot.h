#ifndef PHONONSQUEEZING_SETUPBASICOPERATORSONLYPHOT_CPP_H
#define PHONONSQUEEZING_SETUPBASICOPERATORSONLYPHOT_CPP_H

#include <vector>
#include <complex>

void setupHelectronicOnlyPhot(std::vector<std::complex<double>> &HElectronic);

void setupAOnlyPhot(std::vector<std::complex<double>> &A);

void setupHCouplingOnlyPhot(std::vector<std::complex<double>> &HCoupling);

void setupDOccOnlyPhot(std::vector<std::complex<double>> &DOcc);

#endif //PHONONSQUEEZING_SETUPBASICOPERATORSONLYPHOT_CPP_H
