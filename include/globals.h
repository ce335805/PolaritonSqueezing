#ifndef PHONONSQUEEZING_GLOBALS_H
#define PHONONSQUEEZING_GLOBALS_H

#include <complex>

constexpr std::complex<double> II = std::complex<double>(0., 1.);
constexpr double PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640;

constexpr double tHop(1.);
constexpr double wPt(2.);
constexpr double U(1.);
constexpr double gPh(0.5);
extern double gE;
constexpr ulong LCHAIN(2);
//constexpr double gE(1.);

//const double wPh(2. * std::sqrt(1. - 4. * gPh / 4.) - 0.25);
//const double wPh(1.24);
const double wPh(1.2360679);
//const double wPh(2. * std::sqrt(1. + 4. * gPh / 4. * 0.390434));
//const double wPh(2. * std::sqrt(1. - 4. * gPh / 4. * 0.109566));
//const double wPh(2.02);
//const double wPh(2.02);


constexpr double wP(0.0);
//extern double wP;

constexpr double wDrive(2.);
constexpr double fDrive(3.);
constexpr int timePointsPerDrivingPeriod (80);
constexpr double dt(2. * PI / wDrive / timePointsPerDrivingPeriod);

constexpr ulong dimElectron(4ul);
//constexpr ulong dimElectron(1 << 2 * LCHAIN);
constexpr ulong dimPhonon(1ul);
constexpr ulong dimPhoton(16ul);
constexpr ulong dimHOnlyPhot(dimElectron * dimPhoton);
constexpr ulong dimHTwoPh(dimElectron * dimPhonon * dimPhonon * dimPhoton);
constexpr ulong dimHOnePh(dimElectron * dimPhonon * dimPhoton);

#endif //PHONONSQUEEZING_GLOBALS_H
