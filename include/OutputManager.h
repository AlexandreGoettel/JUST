//! @file	include/OutputManager.h
//! @brief	Declaration of class to process fit results
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef OUTPUTMANAGER_H_
#define OUTPUTMANAGER_H_

//============================================================================
// Standard includes
#include <memory>
#include <vector>
#include <fstream>
// Project includes
#include "FitResults.h"
#include "DataReader.h"
#include "Parser.h"

//============================================================================
// Method definitions

namespace NuFitter {

void ProcessResults(NuFitData*, NuFitPDFs*, const NuFitConfig, NuFitResults);
void ProcessResults(std::vector<NuFitData*>, NuFitPDFs*, const NuFitConfig, std::vector<NuFitResults>);
void fitToFile(std::ofstream&, NuFitResults, const NuFitConfig);

}  // namespace NuFitter


#endif
