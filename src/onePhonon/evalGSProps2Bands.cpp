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

//double wPh;
//double gE;

void evalGSPropsAsOfWPh() {
  
  double gSqrOverOmega = 10.;
  
  const ulong wSteps(50ul);
  std::vector<double> wArr(wSteps, 0.);
  
  for (ulong ind = 0ul; ind < wSteps; ++ind) {
    wArr[ind] = double(ind + 1ul) / 5. + 5.;
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
    //gE = gSqrOverOmega / std::sqrt(wPh);
  
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
  filename = gsPropName2Bands(gSqrOverOmega);
  
  writeStuffToHdf52Bands(wArr,
                         dOcc0Expectation,
                         dOcc1Expectation,
                         dOccUpDnExpectation,
                         dOccSigSigExpectation,
                         n0Expectation,
                         n1Expectation,
                         filename);
}
