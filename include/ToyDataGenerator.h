//! @file	include/ToyDataGenerator.h
//! @brief	Declaration of class to generate toy data
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef TOYDATAGENERATOR_H_
#define TOYDATAGENERATOR_H_

//============================================================================
// Includes
#include <vector>
#include "Parser.h"

//============================================================================
// Class definition

namespace NuFitter {

// Forward declaration
class NuFitConfig;
class NuFitData;

// NuFitToyData *generateToyData(const NuFitConfig *config);
std::vector<NuFitData*> generateToyData(const NuFitConfig config);

}  // namespace NuFitter


#endif
