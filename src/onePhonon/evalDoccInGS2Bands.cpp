#include <iostream>

#include "globals.h"
#include "matrixOperations.h"
#include "include/onePhonon/setupHGlobal2Bands.h"
#include "include/onePhonon/setupBasicOperators2Bands.h"
#include "evalExpectation.h"

#include <chrono>

void evalDoccInGS2Bands(){
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

  setupGlobalH2Bands(H);
  
  std::vector<std::complex<double>> dOcc0;
  std::vector<std::complex<double>> dOcc1;
  setupDOcc0(dOcc0);
  setupDOcc1(dOcc1);
  
  std::vector<std::complex<double>> dOccInterUpDn;
  std::vector<std::complex<double>> dOccInterSigSig;
  setupInterOrbUpDn(dOccInterUpDn);
  setupInterOrbSigSig(dOccInterSigSig);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Setup of H took " << duration.count() << "ms" << '\n';

  std::cout << "Starting diagonalization ..." << '\n';

  start = std::chrono::high_resolution_clock::now();


  std::vector<double> spectrum = diagonalize(H, dimHOnePh, 'V');

  
  
  std::cout << "GS energy = " << spectrum[0] << '\n';
  std::cout << "GS + 1 energy = " << spectrum[1] << '\n';
  std::cout << "GS + 2 energy = " << spectrum[2] << '\n';
  std::cout << "GS + 3 energy = " << spectrum[3] << '\n';
  std::cout << "GS + 4 energy = " << spectrum[4] << '\n';
  std::cout << "GS + 5 energy = " << spectrum[5] << '\n';

  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Diagonalization took " << duration.count() << "ms" << '\n';

  std::cout << '\n';

  std::vector<std::complex<double>> gs (dimHOnePh, std::complex<double> (0., 0.));
  for(ulong ind = 0ul; ind < dimHOnePh; ++ind){
    gs[ind] = H[ind * dimHOnePh];
  }
  
  double dOccExpec0 = evalExpectation(dOcc0, gs, dimHOnePh);
  double dOccExpec1 = evalExpectation(dOcc1, gs, dimHOnePh);
  double dOccInterUpDnExpec = evalExpectation(dOccInterUpDn, gs, dimHOnePh);
  double dOccInterSigSigExpec = evalExpectation(dOccInterSigSig, gs, dimHOnePh);
  
  std::cout << "<dOcc0> = " << dOccExpec0 << '\n';
  std::cout << "<dOcc1> = " << dOccExpec1 << '\n';
  std::cout << "<dOccInterUpDn> = " << dOccInterUpDnExpec << '\n';
  std::cout << "<dOccInterSigSig> = " << dOccInterSigSigExpec << '\n';
  
  std::cout << "<dOccTot> = " << dOccExpec0 + dOccExpec1 + dOccInterUpDnExpec + dOccInterSigSigExpec << '\n';
}
