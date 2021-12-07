#ifndef PHONONSQUEEZING_TIMESTEP_H
#define PHONONSQUEEZING_TIMESTEP_H

void calcTimeStep(const double tPoint,
                  const double pumpPrefac,
                  const std::vector<std::complex<double>> &H,
                  const std::vector<std::complex<double>> &ODrive,
                  std::vector<std::complex<double>> &state,
                  const ulong dimH);

#endif //PHONONSQUEEZING_TIMESTEP_H
