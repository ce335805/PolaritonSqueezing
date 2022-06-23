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

#endif //PHONONSQUEEZING_SETUPBASICOPERATORS2BANDS_H

