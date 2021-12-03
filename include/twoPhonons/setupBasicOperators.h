#ifndef PHONONSQUEEZING_SETUPBASICOPERATORS_H
#define PHONONSQUEEZING_SETUPBASICOPERATORS_H

void setupHelectronic(std::vector<std::complex<double>> &HElectronic);

void setupB1(std::vector<std::complex<double>> &B1);

void setupB2(std::vector<std::complex<double>> &B2);

void setupA(std::vector<std::complex<double>> &A);

void setupHElectronicSmall(std::vector<std::complex<double>> &HElectronicSmall);

void setupDOcc(std::vector<std::complex<double>> &DOcc);

void setupOccSiteI(std::vector<std::complex<double>> &OccSiteI, const ulong site);

#endif //PHONONSQUEEZING_SETUPBASICOPERATORS_H
