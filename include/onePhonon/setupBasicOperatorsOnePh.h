#ifndef PHONONSQUEEZING_SETUPBASICOPERATORSONEPH_H
#define PHONONSQUEEZING_SETUPBASICOPERATORSONEPH_H

void setupHelectronicOnePh(std::vector<std::complex<double>> &HElectronic);

void setupDOccOnePh(std::vector<std::complex<double>> &DOcc);

void setupOccSiteIOnePh(std::vector<std::complex<double>> &OccSiteI, const ulong site);

void setupDOccOnePhNoPHS(std::vector<std::complex<double>> &DOcc);

void setupB(std::vector<std::complex<double>> &B);

void setupAOnePh(std::vector<std::complex<double>> &A);

#endif //PHONONSQUEEZING_SETUPBASICOPERATORSONEPH_H
