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
          ptInd1 * dimHTwoPh * dimPhonon * dimPhonon * 4ul +
          ph2Ind1 * dimHTwoPh * dimPhonon * 4ul +
          ph1Ind1 * dimHTwoPh * 4ul +
          eInd1 * dimHTwoPh +
          ptInd2 * dimPhonon * dimPhonon * 4ul +
          ph2Ind2 * dimPhonon * 4ul +
          ph1Ind2 * 4ul +
          eInd2;

  //return
  //        ptInd1 * dimHTwoPh * dimPhonon * dimPhonon * 4ul +
  //        ph2Ind1 * dimHTwoPh * 4ul +
  //        ph1Ind1 * dimHTwoPh * dimPhonon * 4ul +
  //        eInd1 * dimHTwoPh +
  //        ptInd2 * dimPhonon * dimPhonon * 4ul +
  //        ph2Ind2 * 4ul +
  //        ph1Ind2 * dimPhonon * 4ul +
  //        eInd2;

  //return
  //        ptInd1 * dimHTwoPh +
  //        ph2Ind1 * dimHTwoPh * dimPhonon * dimPhoton +
  //        ph1Ind1 * dimHTwoPh * dimPhoton +
  //        eInd1 * dimHTwoPh * dimPhonon * dimPhonon * dimPhoton +
  //        ptInd2 +
  //        ph2Ind2 * dimPhonon * dimPhoton +
  //        ph1Ind2 * dimPhoton +
  //        eInd2 * dimPhonon * dimPhonon * dimPhoton;

  //ulong index =
  //        ptInd1 * dimHTwoPh * dimPhonon * dimPhonon +
  //        ph2Ind1 * dimHTwoPh +
  //        ph1Ind1 * dimHTwoPh * dimPhonon +
  //        eInd1 * dimHTwoPh * dimPhonon * dimPhonon * dimPhoton +
  //        ptInd2 * dimPhonon * dimPhonon +
  //        ph2Ind2 +
  //        ph1Ind2 * dimPhonon +
  //        eInd2 * dimPhonon * dimPhonon * dimPhoton;
//
  //if(index >= dimHTwoPh * dimHTwoPh){
  //  std::cout << "Index exceeds boundary!!!" << '\n';
  //}
//
  //return index;
}

inline ulong toGlobalMatrixIndexOne(const ulong ptInd1,
                                    const ulong ptInd2,
                                    const ulong phInd1,
                                    const ulong phInd2,
                                    const ulong eInd1,
                                    const ulong eInd2) {
  return
          ptInd1 * dimHOnePh * dimPhonon * 4ul +
          phInd1 * dimHOnePh * 4ul +
          eInd1 * dimHOnePh +
          ptInd2 * dimPhonon * 4ul +
          phInd2 * 4ul +
          eInd2;

  //return
  //        ptInd1 * dimHOnePh +
  //        phInd1 * dimHOnePh * dimPhoton +
  //        eInd1 * dimHOnePh * dimPhoton * dimPhonon +
  //        ptInd2+
  //        phInd2 * dimPhoton+
  //        eInd2 * dimPhoton * dimPhonon;
}

#endif //PHONONSQUEEZING_UTILS_H
