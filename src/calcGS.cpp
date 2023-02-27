#include <iostream>

#include "globals.h"
#include "matrixOperations.h"
#include "include/onePhonon/setUpGlobalHamiltonianOnePh.h"
#include "include/onePhonon/setupBasicOperatorsOnePh.h"
#include "evalExpectation.h"

#include <chrono>

void calcGS(std::vector<std::complex<double>> &gs, std::vector<std::complex<double>> globalH, const ulong dimH){

  std::cout << "Starting diagonalization ..." << '\n';

  auto start = std::chrono::high_resolution_clock::now();

  std::vector<double> spectrum = diagonalize(globalH, dimH, 'V');

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Diagonalization took " << duration.count() << "ms" << '\n';

  std::cout << '\n';

  for(ulong ind = 0ul; ind < dimH; ++ind){
    gs[ind] = globalH[ind * dimH];
  }

}

double calcGSWithE(std::vector<std::complex<double>> &gs, std::vector<std::complex<double>> globalH, const ulong dimH){
  
  std::cout << "Starting diagonalization ..." << '\n';
  
  auto start = std::chrono::high_resolution_clock::now();
  
  std::vector<double> spectrum = diagonalize(globalH, dimH, 'V');
  
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  
  std::cout << "Diagonalization took " << duration.count() << "ms" << '\n';
  
  std::cout << '\n';
  
  for(ulong ind = 0ul; ind < dimH; ++ind){
    gs[ind] = globalH[ind * dimH];
  }
  
  return spectrum[0];
  
}

std::vector<double> calcEigenEnergies(std::vector<std::complex<double>> &globalH, const ulong dimH){
  
  std::cout << "Starting diagonalization ..." << '\n';
  
  auto start = std::chrono::high_resolution_clock::now();
  std::vector<double> spectrum = diagonalize(globalH, dimH, 'V');
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  
  std::cout << "Diagonalization took " << duration.count() << "ms" << '\n';
  std::cout << '\n';
  
  return spectrum;
  
}
