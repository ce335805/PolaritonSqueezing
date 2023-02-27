#ifndef PHONONSQUEEZING_SETUPBASICOPERATORS2BANDS_H
#define PHONONSQUEEZING_SETUPBASICOPERATORS2BANDS_H

#include "vector"
#include "complex"

void setupHelectronic2Bands(std::vector<std::complex<double>> &HElectronic);

void setupDOcc0(std::vector<std::complex<double>> &DOcc);

void setupDOcc1(std::vector<std::complex<double>> &DOcc);

void setupInterbandHop0(std::vector<std::complex<double>> &hopInter);

void setupInterbandHop1(std::vector<std::complex<double>> &hopInter);

void setupA0(std::vector<std::complex<double>> &A0);

void setupA1(std::vector<std::complex<double>> &A1);

void setupInterOrbUpDn(std::vector<std::complex<double>> &interOrbUpDn);

void setupInterOrbSigSig(std::vector<std::complex<double>> &interOrbSigSig);

void setupN0(std::vector<std::complex<double>> &N0);

void setupN1(std::vector<std::complex<double>> &N1);

void setupNC0(std::vector<std::complex<double>> &NC0);
void setupND0(std::vector<std::complex<double>> &ND0);
void setupNC1(std::vector<std::complex<double>> &NC1);
void setupND1(std::vector<std::complex<double>> &ND1);

#endif //PHONONSQUEEZING_SETUPBASICOPERATORS2BANDS_H

