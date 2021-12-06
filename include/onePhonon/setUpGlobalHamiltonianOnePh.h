#ifndef PHONONSQUEEZING_SETUPGLOBALHAMILTONIANONEPH_H
#define PHONONSQUEEZING_SETUPGLOBALHAMILTONIANONEPH_H

#include <vector>
#include <complex>

#include "globals.h"

void setupGlobalHOnePh(std::vector<std::complex<double>> &H);

void setupPhPtCouplingOnePh(std::vector<std::complex<double>> &HPhPt);

void setupQuadEPhCouplingOnePh(std::vector<std::complex<double>> &HEPh);

void setupUncoupledHamiltonianOnePh(std::vector<std::complex<double>> &H);

#endif //PHONONSQUEEZING_SETUPGLOBALHAMILTONIANONEPH_H
