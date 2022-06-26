#ifndef PHONONSQUEEZING_SETUPHGLOBAL2BANDS_H
#define PHONONSQUEEZING_SETUPHGLOBAL2BANDS_H

#include <vector>
#include <complex>

#include "globals.h"

void setupGlobalH2Bands(std::vector<std::complex<double>> &H);
void setupHEBos(std::vector<std::complex<double>> &HEBos);
void setupHFreePhon(std::vector<std::complex<double>> &freePhonon);

#endif //PHONONSQUEEZING_SETUPHGLOBAL2BANDS_H
