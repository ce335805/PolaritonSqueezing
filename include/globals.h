#ifndef PHONONSQUEEZING_GLOBALS_H
#define PHONONSQUEEZING_GLOBALS_H

#include <complex>

constexpr std::complex<double> II = std::complex<double>(0., 1.);
constexpr double PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640;

constexpr double tHop(0.25);
constexpr double wPt(0.1);
constexpr double wPh(0.2);
constexpr double U(5. * tHop);
constexpr double gPh(-0.1);
constexpr double wP(0.2);

constexpr double wDrive(0.15);
constexpr double fDrive(.5);
constexpr double dt(2. * PI / wDrive / 40.);

constexpr ulong dimPhonon(14ul);
constexpr ulong dimPhoton(14ul);
constexpr ulong dimHTwoPh(4ul * dimPhonon * dimPhonon * dimPhoton);
constexpr ulong dimHOnePh(4ul * dimPhonon * dimPhoton);


#endif //PHONONSQUEEZING_GLOBALS_H
