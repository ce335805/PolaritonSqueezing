#ifndef PHONONSQUEEZING_SETUPELECTRONICOPERATORSSMALL_CPP_H
#define PHONONSQUEEZING_SETUPELECTRONICOPERATORSSMALL_CPP_H

void setupHElectronicSmall(std::vector<std::complex<double>> &HElectronicSmall);

void setupDoubleOccSmall(std::vector<std::complex<double>> &DOccSmall);

void setupDoubleOccSmallNoPHS(std::vector<std::complex<double>> &DOccSmall);

void setupDOccSiteISmall(std::vector<std::complex<double>> &dOccSmallSiteI, const ulong site);

void setupOccSiteISmall(std::vector<std::complex<double>> &occSmallSiteI, const ulong site);

#endif //PHONONSQUEEZING_SETUPELECTRONICOPERATORSSMALL_CPP_H
