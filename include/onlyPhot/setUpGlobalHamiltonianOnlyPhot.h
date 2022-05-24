#ifndef PHONONSQUEEZING_SETUPGLOBALHAMILTONIANONLYPHOT_H
#define PHONONSQUEEZING_SETUPGLOBALHAMILTONIANONLYPHOT_H

#include <vector>
#include <complex>

void setupGlobalHOnlyPhot(std::vector<std::complex<double>> &H);

void setupUncoupledHamiltonianOnlyPhot(std::vector<std::complex<double>> &H);

void setupEPtCoupling(std::vector<std::complex<double>> &HEPt);

#endif //PHONONSQUEEZING_SETUPGLOBALHAMILTONIANONLYPHOT_H
