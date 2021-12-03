#include <iostream>

#include "globals.h"
#include "matrixOperations.h"
#include "include/twoPhonons/setUpGlobalHamiltonian.cpp.h"
#include "include/twoPhonons/setupBasicOperators.h"
#include "include/twoPhonons/evalExpectation.h"

#include <chrono>


int main() {

  std::cout << "Hello World! - let's squeeze some phonons" << std::endl;

  std::cout << "Using parameters: " << '\n';
  std::cout << "w-Phonon = " << wPh << '\n';
  std::cout << "w-Photon = " << wPt << '\n';
  std::cout << "electronic U = " << U << '\n';
  std::cout << "gPh = " << gPh << '\n';
  std::cout << "wP = " << wP << '\n';

  std::cout << "\n";

  std::vector<std::complex<double>> H;
  //setupUncoupledHamiltonian(H);

  std::cout << "Starting setup of Hamiltonian ..." << '\n';

  auto start = std::chrono::high_resolution_clock::now();

  setupGlobalH(H);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  auto timeForDiagonalization = duration.count();

  std::cout << "Setup of H took " << duration.count() << "ms" << '\n';

  //for(ulong ind1 = 0ul; ind1 < dimHOnePh; ++ind1){
  //  for(ulong ind2 = 0ul; ind2 < dimHOnePh; ++ind2) {
  //  std::cout << H[ind1 * dimHOnePh + ind2] << " ";
  //  }
  //  std::cout << '\n';
  //}

  std::cout << "Starting diagonalization ..." << '\n';

  start = std::chrono::high_resolution_clock::now();

  std::vector<double> spectrum = diagonalize(H, dimHOnePh, 'V');

  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  timeForDiagonalization = duration.count();

  std::cout << "Diagonalization took " << duration.count() << "ms" << '\n';

  std::cout << '\n';

  //for(ulong ind = 0ul; ind < dimHOnePh; ++ind){
  //  std::cout << spectrum[ind] << '\n';
  //}

  std::cout << '\n';

  std::vector<std::complex<double>> gs (dimHOnePh, std::complex<double> (0., 0.));
  for(ulong ind = 0ul; ind < dimHOnePh; ++ind){
    gs[ind] = H[ind * dimHOnePh];
  }

  //for(ulong ind = 0ul; ind < dimHOnePh; ++ind){
  //  std::cout << gs[ind] << '\n';
  //}

  std::vector<std::complex<double>> dOcc;
  setupDOcc(dOcc);

  double dOccExpec = evalExpectation(dOcc, gs);

  std::cout << "<dOcc> = " << dOccExpec + 0.5 << '\n';
  std::cout << "Delta<dOcc> = " << (dOccExpec + 0.5) - 0.109566 << '\n';

  return 0;
}

