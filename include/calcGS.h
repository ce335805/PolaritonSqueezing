#ifndef PHONONSQUEEZING_CALCGS_H
#define PHONONSQUEEZING_CALCGS_H

void calcGS(std::vector<std::complex<double>> &gs, std::vector<std::complex<double>> globalH, const ulong dimH);

double calcGSWithE(std::vector<std::complex<double>> &gs, std::vector<std::complex<double>> globalH, const ulong dimH);

std::vector<double> calcEigenEnergies(std::vector<std::complex<double>> &globalH, const ulong dimH);

#endif //PHONONSQUEEZING_CALCGS_H
