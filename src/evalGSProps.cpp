#include <vector>
#include <complex>
#include <iostream>


#include "globals.h"
#include "setupOps.h"
#include "evalExpectation.h"
#include "setUpGlobalHamiltonianOnePh.h"
#include "setUpGlobalHamiltonian.h"
#include "calcGS.h"
#include "writeStuffToHdf5.h"
#include "makeFilenames.h"

double wP;

void evalGSProps(const bool twoPhonons){

  const ulong dimH = twoPhonons ? dimHTwoPh : dimHOnePh;

  const ulong wPSteps (11ul);
  std::vector<double> wPArr (wPSteps, 0.);

  for(ulong ind = 0ul; ind < wPSteps; ++ind){
    wPArr[ind] = double(ind) / 2.;
  }

  std::vector<std::complex<double>> dOcc;
  std::vector<std::complex<double>> Xpt;
  std::vector<std::complex<double>> XptSqr;
  std::vector<std::complex<double>> Npt;
  std::vector<std::complex<double>> X1ph;
  std::vector<std::complex<double>> X1phSqr;
  std::vector<std::complex<double>> N1ph;
  std::vector<std::complex<double>> X2ph;
  std::vector<std::complex<double>> X2phSqr;
  std::vector<std::complex<double>> N2ph;

  if (twoPhonons) {
    setupOpsTwoPh(dOcc, Xpt, XptSqr, Npt, X1ph, X1phSqr, N1ph, X2ph, X2phSqr, N2ph);
  } else {
    setupOpsOnePh(dOcc, Xpt, XptSqr, Npt, X1ph, X1phSqr, N1ph, X2ph, X2phSqr, N2ph);
  }

  std::vector<double> dOccExpectation(wPSteps, 0.);
  std::vector<double> XptExpectation(wPSteps, 0.);
  std::vector<double> XptSqrExpectation(wPSteps, 0.);
  std::vector<double> NptExpectation(wPSteps, 0.);
  std::vector<double> X1phExpectation(wPSteps, 0.);
  std::vector<double> X1phSqrExpectation(wPSteps, 0.);
  std::vector<double> N1phExpectation(wPSteps, 0.);
  std::vector<double> X2phExpectation;
  std::vector<double> X2phSqrExpectation;
  std::vector<double> N2phExpectation;
  if (twoPhonons) {
    X2phExpectation = std::vector<double>(wPSteps, 0.);
    X2phSqrExpectation = std::vector<double>(wPSteps, 0.);
    N2phExpectation = std::vector<double>(wPSteps, 0.);
  }

  std::vector<std::complex<double>> gs(dimH, std::complex<double> (0., 0.));
  std::vector<std::complex<double>> H;

  for (ulong wPStep = 0ul; wPStep < wPSteps; ++wPStep) {

    wP = wPArr[wPStep];
    std::cout << "wP = " << wP << '\n';
    if(twoPhonons){
      setupGlobalH(H);
    } else {
      setupGlobalHOnePh(H);
    }

    calcGS(gs, H, dimH);

    dOccExpectation[wPStep] = evalExpectation(dOcc, gs, dimH);
    XptExpectation[wPStep] = evalExpectation(Xpt, gs, dimH);
    XptSqrExpectation[wPStep] = evalExpectation(XptSqr, gs, dimH);
    NptExpectation[wPStep] = evalExpectation(Npt, gs, dimH);
    X1phExpectation[wPStep] = evalExpectation(X1ph, gs, dimH);
    X1phSqrExpectation[wPStep] = evalExpectation(X1phSqr, gs, dimH);
    N1phExpectation[wPStep] = evalExpectation(N1ph, gs, dimH);
    if (twoPhonons) {
      X2phExpectation[wPStep] = evalExpectation(X2ph, gs, dimH);
      X2phSqrExpectation[wPStep] = evalExpectation(X2phSqr, gs, dimH);
      N2phExpectation[wPStep] = evalExpectation(N2ph, gs, dimH);
    }
  }

  std::string filename;
  filename = gsPropName(twoPhonons);

  writeStuffToHdf5(wPArr,
                   wPArr,
                   dOccExpectation,
                   XptExpectation,
                   XptSqrExpectation,
                   NptExpectation,
                   X1phExpectation,
                   X1phSqrExpectation,
                   N1phExpectation,
                   X2phExpectation,
                   X2phSqrExpectation,
                   N2phExpectation,
                   filename,
                   twoPhonons);


}

