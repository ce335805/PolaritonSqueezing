#ifndef PHONONSQUEEZING_UTILS_H
#define PHONONSQUEEZING_UTILS_H

#include "globals.h"


inline ulong toGlobalMatrixIndex(const ulong ptInd1,
                                 const ulong ptInd2,
                                 const ulong ph2Ind1,
                                 const ulong ph2Ind2,
                                 const ulong ph1Ind1,
                                 const ulong ph1Ind2,
                                 const ulong eInd1,
                                 const ulong eInd2) {
  return
          ptInd1 * dimHTwoPh * dimPhonon * dimPhonon * dimElectron +
          ph2Ind1 * dimHTwoPh * dimPhonon * dimElectron +
          ph1Ind1 * dimHTwoPh * dimElectron +
          eInd1 * dimHTwoPh +
          ptInd2 * dimPhonon * dimPhonon * dimElectron +
          ph2Ind2 * dimPhonon * dimElectron +
          ph1Ind2 * dimElectron +
          eInd2;
}

inline ulong toGlobalMatrixIndexOne(const ulong ptInd1,
                                    const ulong ptInd2,
                                    const ulong phInd1,
                                    const ulong phInd2,
                                    const ulong eInd1,
                                    const ulong eInd2) {
  return
          ptInd1 * dimHOnePh * dimPhonon * dimElectron +
          phInd1 * dimHOnePh * dimElectron +
          eInd1 * dimHOnePh +
          ptInd2 * dimPhonon * dimElectron +
          phInd2 * dimElectron +
          eInd2;
}

inline ulong toGlobalMatrixIndexOnlyPhot(const ulong ptInd1,
                                    const ulong ptInd2,
                                    const ulong eInd1,
                                    const ulong eInd2) {
  return
      ptInd1 * dimHOnlyPhot * dimElectron +
      eInd1 * dimHOnlyPhot +
      ptInd2 * dimElectron +
      eInd2;
}

#endif //PHONONSQUEEZING_UTILS_H
