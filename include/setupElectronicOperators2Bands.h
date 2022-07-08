#ifndef PHONONSQUEEZING_SETUPELECTRONICOPERATORS2BANDS_H
#define PHONONSQUEEZING_SETUPELECTRONICOPERATORS2BANDS_H

#include <vector>
#include <complex>

void setupHElectronicSmall2Bands(std::vector<std::complex<double>> &HElectronicSmall);

void setupDoubleOccSmall0(std::vector<std::complex<double>> &DOccSmall);

void setupDoubleOccSmall1(std::vector<std::complex<double>> &DOccSmall);

void setupInterbandHop0Small(std::vector<std::complex<double>> &interbandHop);

void setupInterbandHop1Small(std::vector<std::complex<double>> &interbandHop);

void setupInterbandHop1Small(std::vector<std::complex<double>> &interbandHop);

void setupInterOrbUpDnSmall(std::vector<std::complex<double>> &interOrbUpDn);

void setupInterOrbSigSigSmall(std::vector<std::complex<double>> &interOrbSigSig);

void setupN0Small(std::vector<std::complex<double>> &N0);

void setupN1Small(std::vector<std::complex<double>> &N1);

#endif //PHONONSQUEEZING_SETUPELECTRONICOPERATORS2BANDS_H
