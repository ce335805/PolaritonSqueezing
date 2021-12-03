#include <iostream>

#include "globals.h"
#include "matrixOperations.h"
#include "include/onePhonon/setUpGlobalHamiltonianOnePh.h"
#include "include/onePhonon/setupBasicOperatorsOnePh.h"
#include "include/onePhonon/evalExpectationOnePh.h"

#include <chrono>

void evalDoccInGSOnePh(){
  std::cout << "Using parameters: " << '\n';
  std::cout << "w-Phonon = " << wPh << '\n';
  std::cout << "w-Photon = " << wPt << '\n';
  std::cout << "electronic U = " << U << '\n';
  std::cout << "gPh = " << gPh << '\n';
  std::cout << "wP = " << wP << '\n';

  std::cout << "\n";

  std::vector<std::complex<double>> H;

  std::cout << "Starting setup of Hamiltonian ..." << '\n';

  auto start = std::chrono::high_resolution_clock::now();

  setupGlobalHOnePh(H);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Setup of H took " << duration.count() << "ms" << '\n';

  std::cout << "Starting diagonalization ..." << '\n';

  start = std::chrono::high_resolution_clock::now();

  std::vector<double> spectrum = diagonalize(H, dimHOnePh, 'V');

  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Diagonalization took " << duration.count() << "ms" << '\n';

  std::cout << '\n';

  std::vector<std::complex<double>> gs (dimHOnePh, std::complex<double> (0., 0.));
  for(ulong ind = 0ul; ind < dimHOnePh; ++ind){
    gs[ind] = H[ind * dimHOnePh];
  }

  std::vector<std::complex<double>> dOcc;
  setupDOccOnePh(dOcc);

  double dOccExpec = evalExpectationOnePh(dOcc, gs);

  std::cout << "<dOcc> = " << dOccExpec + 0.5 << '\n';
  std::cout << "Delta<dOcc> = " << (dOccExpec + 0.5) - 0.109566 << '\n';
}
