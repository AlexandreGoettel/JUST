//! @file	include/Fitter.h
//! @brief	Declaration of class to fit MC PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef FITTER_H_
#define FITTER_H_

//============================================================================
// Includes
#include "FitResults.h"



//============================================================================
// Method definitions

namespace NuFitter {

//============================================================================
// Forward declarations
class NuFitConfig;
class NuFitPDFs;
class NuFitData;
class NuFitToyData;

namespace MCFit {

// Overloaded functions to fit Data or Toy Data
NuFitResults *Fit(NuFitData*, NuFitPDFs*, const NuFitConfig*);
NuFitResults *Fit(NuFitToyData*, NuFitPDFs*, const NuFitConfig*);

}  // namespace MCFit
}  // namespace NuFitter


#endif
