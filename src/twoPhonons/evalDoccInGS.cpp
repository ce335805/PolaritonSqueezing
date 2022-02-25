#include <iostream>

#include "globals.h"
#include "matrixOperations.h"
#include "include/twoPhonons/setUpGlobalHamiltonian.h"
#include "include/twoPhonons/setupBasicOperators.h"
#include "include/evalExpectation.h"

#include <chrono>

void evalDoccInGSTwoPh(){
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

  setupGlobalH(H);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Setup of H took " << duration.count() << "ms" << '\n';

  std::cout << "Starting diagonalization ..." << '\n';

  start = std::chrono::high_resolution_clock::now();

  std::vector<double> spectrum = diagonalize(H, dimHTwoPh, 'V');

  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Diagonalization took " << duration.count() << "ms" << '\n';

  std::cout << '\n';

  std::vector<std::complex<double>> gs (dimHTwoPh, std::complex<double> (0., 0.));
  for(ulong ind = 0ul; ind < dimHTwoPh; ++ind){
    gs[ind] = H[ind * dimHTwoPh];
  }

  std::vector<std::complex<double>> dOcc;
  setupDOcc(dOcc);

  double dOccExpec = evalExpectation(dOcc, gs, dimHTwoPh);

  std::cout << "<dOcc> = " << dOccExpec + 0.5 << '\n';
  std::cout << "Delta<dOcc> = " << (dOccExpec + 0.5) - 0.109566 << '\n';
}

void evalDoccThermal(const double temp){
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

  setupGlobalH(H);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Setup of H took " << duration.count() << "ms" << '\n';

  std::cout << "Starting diagonalization ..." << '\n';

  start = std::chrono::high_resolution_clock::now();

  std::vector<double> spectrum = diagonalize(H, dimHTwoPh, 'V');

  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Diagonalization took " << duration.count() << "ms" << '\n';

  std::cout << '\n';

  std::vector<std::complex<double>> gs (dimHTwoPh, std::complex<double> (0., 0.));
  for(ulong ind = 0ul; ind < dimHTwoPh; ++ind){
    gs[ind] = H[ind * dimHTwoPh];
  }

  std::vector<std::complex<double>> dOcc;
  setupDOcc(dOcc);

  double dOccExpec = evalExpectation(dOcc, gs, dimHTwoPh);

  std::cout << "<dOcc> = " << dOccExpec + 0.5 << '\n';
  std::cout << "Delta<dOcc> = " << (dOccExpec + 0.5) - 0.109566 << '\n';
}
