//! @file	include/OutputManager.h
//! @brief	Declaration of class to process fit results
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef OUTPUTMANAGER_H_
#define OUTPUTMANAGER_H_

//============================================================================
// Includes

//============================================================================
// Method definitions

namespace NuFitter {

//============================================================================
// Forward declarations
class NuFitResults;

void ProcessResults(NuFitResults*);

}  // namespace NuFitter


#endif
