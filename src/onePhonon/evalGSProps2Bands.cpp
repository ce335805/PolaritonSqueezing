#include <iostream>

#include "globals.h"
#include "matrixOperations.h"
#include "evalExpectation.h"
#include "setupHGlobal2Bands.h"
#include "setupBasicOperators2Bands.h"
#include "setupOps.h"
#include "calcGS.h"
#include "writeStuffToHdf5.h"
#include "makeFilenames.h"

#include <chrono>

double wPh;
//double gE;

void evalGSPropsAsOfWPh() {
  
  const ulong wSteps(100ul);
  std::vector<double> wArr(wSteps, 0.);
  
  for (ulong ind = 0ul; ind < wSteps; ++ind) {
    wArr[ind] = double(ind + 1ul) / 10. + 5.;
  }
  
  std::vector<std::complex<double>> dOcc0;
  std::vector<std::complex<double>> dOcc1;
  std::vector<std::complex<double>> dOccUpDn;
  std::vector<std::complex<double>> dOccSigSig;
  std::vector<std::complex<double>> n0;
  std::vector<std::complex<double>> n1;
  
  setupOps2Bands(dOcc0, dOcc1, dOccUpDn, dOccSigSig, n0, n1);
  
  std::vector<double> dOcc0Expectation(wSteps, 0.);
  std::vector<double> dOcc1Expectation(wSteps, 0.);
  std::vector<double> dOccUpDnExpectation(wSteps, 0.);
  std::vector<double> dOccSigSigExpectation(wSteps, 0.);
  std::vector<double> eGSExpectation(wSteps, 0.);
  std::vector<double> n0Expectation(wSteps, 0.);
  std::vector<double> n1Expectation(wSteps, 0.);
  
  std::vector<std::complex<double>> gs(dimHOnePh, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> H;
  
  for (ulong wStep = 0ul; wStep < wSteps; ++wStep) {
    
    ////////////////// set wPh //////////////////////
    //wPh = wArr[wStep];
  
    std::cout << "wPh = " << wPh << '\n';
    std::cout << "gE = " << gE << '\n';
    setupGlobalH2Bands(H);
  
    //eGSExpectation[wStep] = calcGSWithE(gs, H, dimHOnePh);
    calcGS(gs, H, dimHOnePh);
    dOcc0Expectation[wStep] = evalExpectation(dOcc0, gs, dimHOnePh);
    dOcc1Expectation[wStep] = evalExpectation(dOcc1, gs, dimHOnePh);
    dOccUpDnExpectation[wStep] = evalExpectation(dOccUpDn, gs, dimHOnePh);
    dOccSigSigExpectation[wStep] = evalExpectation(dOccSigSig, gs, dimHOnePh);
    n0Expectation[wStep] = evalExpectation(n0, gs, dimHOnePh);
    n1Expectation[wStep] = evalExpectation(n1, gs, dimHOnePh);
  }
  
  std::string filename;
  filename = gsPropName2Bands(gE);
  
  writeStuffToHdf52Bands(wArr,
                         dOcc0Expectation,
                         dOcc1Expectation,
                         dOccUpDnExpectation,
                         dOccSigSigExpectation,
                         n0Expectation,
                         n1Expectation,
                         filename);
}


void eigenenergiesAsOfG() {
  
  const ulong gSteps(300ul);
  std::vector<double> gArr(gSteps, 0.);
  for (ulong ind = 0ul; ind < gSteps; ++ind) {
    gArr[ind] = double(ind) / 100.;
  }
  
  std::vector<std::complex<double>> nc0;
  std::vector<std::complex<double>> nd0;
  std::vector<std::complex<double>> nc1;
  std::vector<std::complex<double>> nd1;
  std::vector<std::complex<double>> nBos0;
  std::vector<std::complex<double>> nBos1;
  
  setupOps2BandsSpectrum(nc0, nd0, nc1, nd1, nBos0, nBos1);
  
  std::vector<double> nc0Expectation(gSteps * dimHOnePh, 0.);
  std::vector<double> nd0Expectation(gSteps * dimHOnePh, 0.);
  std::vector<double> nc1Expectation(gSteps * dimHOnePh, 0.);
  std::vector<double> nd1Expectation(gSteps * dimHOnePh, 0.);
  std::vector<double> nBos0Expectation(gSteps * dimHOnePh,0.);
  std::vector<double> nBos1Expectation(gSteps * dimHOnePh,0.);
  
  std::vector<double> nc0ExpectationTemp(dimHOnePh, 0.);
  std::vector<double> nd0ExpectationTemp(dimHOnePh, 0.);
  std::vector<double> nc1ExpectationTemp(dimHOnePh, 0.);
  std::vector<double> nd1ExpectationTemp(dimHOnePh, 0.);
  std::vector<double> nBos0ExpectationTemp(dimHOnePh, 0.);
  std::vector<double> nBos1ExpectationTemp(dimHOnePh, 0.);
  
  std::vector<std::complex<double>> gs(dimHOnePh, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> ex(dimHOnePh, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> H;
  std::vector<double> spectrum(gSteps * dimHOnePh, 0.);
  
  setupGlobalH2Bands(H);
  
  for (ulong gStep = 0ul; gStep < gSteps; ++gStep) {
    
    ////////////////// set gE //////////////////////
    //gE = gArr[gStep];
    std::cout << "gE = " << gE << '\n';
    setupGlobalH2Bands(H);
    std::vector<double> spectrumTemp = calcEigenEnergies(H, dimHOnePh);
  
    for(ulong ind = 0ul; ind < dimHOnePh; ++ind){
      gs[ind] = H[ind * dimHOnePh];
      ex[ind] = H[ind * dimHOnePh + 1];
    }
    
    //std::cout << spectrumTemp[0] << '\n';
    //std::cout << spectrumTemp[1] << '\n';
  
  
    double nBos0_gs = evalExpectation(nBos0, gs, dimHOnePh);
    double nBos0_ex = evalExpectation(nBos0, ex, dimHOnePh);
    
    double nBos1_gs = evalExpectation(nBos1, gs, dimHOnePh);
    double nBos1_ex = evalExpectation(nBos1, ex, dimHOnePh);
  
    //std::cout << nBos0_gs + nBos0_ex << '\n';
    //std::cout << nBos1_gs + nBos1_ex << '\n';
    
    evalAllExpectations(nc0ExpectationTemp, nc0, H, dimHOnePh);
    evalAllExpectations(nd0ExpectationTemp, nd0, H, dimHOnePh);
    evalAllExpectations(nc1ExpectationTemp, nc1, H, dimHOnePh);
    evalAllExpectations(nd1ExpectationTemp, nd1, H, dimHOnePh);
    evalAllExpectations(nBos0ExpectationTemp, nBos0, H, dimHOnePh);
    evalAllExpectations(nBos1ExpectationTemp, nBos1, H, dimHOnePh);
  
    std::cout << nBos0ExpectationTemp[0] + nBos0ExpectationTemp[1] << '\n';
    std::cout << nBos1ExpectationTemp[0] + nBos1ExpectationTemp[1] << '\n';
    
    for(ulong hInd = 0ul; hInd < dimHOnePh; ++hInd){
      spectrum[dimHOnePh * gStep + hInd] = spectrumTemp[hInd];
      nc0Expectation[dimHOnePh * gStep + hInd] = nc0ExpectationTemp[hInd];
      nd0Expectation[dimHOnePh * gStep + hInd] = nd0ExpectationTemp[hInd];
      nc1Expectation[dimHOnePh * gStep + hInd] = nc1ExpectationTemp[hInd];
      nd1Expectation[dimHOnePh * gStep + hInd] = nd1ExpectationTemp[hInd];
      nBos0Expectation[dimHOnePh * gStep + hInd] = nBos0ExpectationTemp[hInd];
      nBos1Expectation[dimHOnePh * gStep + hInd] = nBos1ExpectationTemp[hInd];
    }
    //std::cout << spectrum[dimHOnePh * gStep] << '\n';
    //std::cout << nBos0Expectation[dimHOnePh * gStep] + nBos0Expectation[dimHOnePh * gStep + 1] << '\n';
    //std::cout << nBos1Expectation[dimHOnePh * gStep] + nBos1Expectation[dimHOnePh * gStep + 1] << '\n';
  }
  
  std::string filename;
  filename = spectrumName();
  
  writeSpectrumToFile(gArr,
                      spectrum,
                      nc0Expectation,
                      nd0Expectation,
                      nc1Expectation,
                      nd1Expectation,
                      nBos0Expectation,
                      nBos1Expectation,
                      filename);
}
