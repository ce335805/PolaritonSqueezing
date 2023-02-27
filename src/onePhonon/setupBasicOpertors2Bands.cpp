#include <complex>
#include <vector>
#include <iostream>

#include "globals.h"
#include "utils.h"
#include "setupElectronicOperators2Bands.h"

void setupHelectronic2Bands(std::vector<std::complex<double>> &HElectronic) {

  HElectronic = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> HElectronicSmall;
  setupHElectronicSmall2Bands(HElectronicSmall);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong ph1Ind = 0ul; ph1Ind < dimPhonon; ++ph1Ind) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          HElectronic[toGlobalMatrixIndexOne(ptInd, ptInd, ph1Ind, ph1Ind, eInd1, eInd2)] = HElectronicSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }

}

void setupDOcc0(std::vector<std::complex<double>> &DOcc) {

  DOcc = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> DOccSmall;
  setupDoubleOccSmall0(DOccSmall);

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          DOcc[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)] = DOccSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupDOcc1(std::vector<std::complex<double>> &DOcc) {
  
  DOcc = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> DOccSmall;
  setupDoubleOccSmall1(DOccSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          DOcc[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)] = DOccSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupInterbandHop0(std::vector<std::complex<double>> &hopInter) {
  
  hopInter = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> hopInterSmall;
  setupInterbandHop0Small(hopInterSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          hopInter[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
          = hopInterSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupInterbandHop1(std::vector<std::complex<double>> &hopInter) {
  
  hopInter = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> hopInterSmall;
  setupInterbandHop1Small(hopInterSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          hopInter[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
              = hopInterSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupA0(std::vector<std::complex<double>> &A0) {
  A0 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd1 = 0ul; phInd1 < dimPhonon; ++phInd1) {
      for (ulong phInd2 = 0ul; phInd2 < dimPhonon; ++phInd2) {
        for (ulong eInd = 0ul; eInd < dimElectron; ++eInd) {
          if (phInd2 == phInd1 + 1ul) {
            A0[toGlobalMatrixIndexOne(ptInd, ptInd, phInd1, phInd2, eInd, eInd)] = std::sqrt(double(phInd2));
          }
        }
      }
    }
  }

}

void setupA1(std::vector<std::complex<double>> &A1) {
  A1 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  for (ulong ptInd1 = 0ul; ptInd1 < dimPhoton; ++ptInd1) {
    for (ulong ptInd2 = 0ul; ptInd2 < dimPhoton; ++ptInd2) {
      for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
        for (ulong eInd = 0ul; eInd < dimElectron; ++eInd) {
          if (ptInd2 == ptInd1 + 1ul) {
            A1[toGlobalMatrixIndexOne(ptInd1, ptInd2, phInd, phInd, eInd, eInd)] = std::sqrt(double(ptInd2));
          }
        }
      }
    }
  }

}

void setupInterOrbUpDn(std::vector<std::complex<double>> &interOrbUpDn) {
  
  interOrbUpDn = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> interOrbUpDnSmall;
  setupInterOrbUpDnSmall(interOrbUpDnSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          interOrbUpDn[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
              = interOrbUpDnSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupInterOrbSigSig(std::vector<std::complex<double>> &interOrbSigSig) {
  
  interOrbSigSig = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> interOrbSigSigSmall;
  setupInterOrbSigSigSmall(interOrbSigSigSmall);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          interOrbSigSig[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
              = interOrbSigSigSmall[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupN0(std::vector<std::complex<double>> &N0) {
  
  N0 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> n0small;
  setupN0Small(n0small);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          N0[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
              = n0small[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}


void setupN1(std::vector<std::complex<double>> &N1) {
  
  N1 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> n1small;
  setupNC0Small(n1small);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          N1[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
              = n1small[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupNC0(std::vector<std::complex<double>> &NC0) {
  
  NC0 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> nc0small;
  setupNC0Small(nc0small);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          NC0[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
              = nc0small[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupND0(std::vector<std::complex<double>> &ND0) {
  
  ND0 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> nd0small;
  setupND0Small(nd0small);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          ND0[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
              = nd0small[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupNC1(std::vector<std::complex<double>> &NC1) {
  
  NC1 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> nc1small;
  setupNC1Small(nc1small);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          NC1[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
              = nc1small[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}

void setupND1(std::vector<std::complex<double>> &ND1) {
  
  ND1 = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> nd1small;
  setupND1Small(nd1small);
  
  for (ulong ptInd = 0ul; ptInd < dimPhoton; ++ptInd) {
    for (ulong phInd = 0ul; phInd < dimPhonon; ++phInd) {
      for (ulong eInd1 = 0ul; eInd1 < dimElectron; ++eInd1) {
        for (ulong eInd2 = 0ul; eInd2 < dimElectron; ++eInd2) {
          ND1[toGlobalMatrixIndexOne(ptInd, ptInd, phInd, phInd, eInd1, eInd2)]
              = nd1small[eInd1 * dimElectron + eInd2];
        }
      }
    }
  }
}