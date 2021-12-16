#ifndef PHONONSQUEEZING_GLOBALS_H
#define PHONONSQUEEZING_GLOBALS_H

#include <complex>

constexpr std::complex<double> II = std::complex<double>(0., 1.);
constexpr double PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640;

constexpr double tHop(1.);
constexpr double wPt(2.);
constexpr double U(5.);
constexpr double gPh(0.5);

const double wPh(2. * std::sqrt(1. - 4. * gPh / 4.));
//const double wPh(2. * std::sqrt(1. + 4. * gPh / 4. * 0.390434));
//const double wPh(2. * std::sqrt(1. - 4. * gPh / 4. * 0.109566));


//constexpr double wP(0.5);
extern double wP;

constexpr double wDrive(2.);
constexpr double fDrive(2.);
constexpr int timePointsPerDrivingPeriod (80);
constexpr double dt(2. * PI / wDrive / timePointsPerDrivingPeriod);

constexpr ulong dimPhonon(6ul);
constexpr ulong dimPhoton(6ul);
constexpr ulong dimHTwoPh(4ul * dimPhonon * dimPhonon * dimPhoton);
constexpr ulong dimHOnePh(4ul * dimPhonon * dimPhoton);

#endif //PHONONSQUEEZING_GLOBALS_H
