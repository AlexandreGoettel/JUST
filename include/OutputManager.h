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
// ROOT includes
#include "TFile.h"
// Project includes
#include "FitResults.h"
#include "DataReader.h"
#include "Parser.h"

//============================================================================
// Method definitions

namespace NuFitter {

// TODO: maybe a class here would be useful after all?
void ProcessResults(NuFitData*, NuFitPDFs*, const NuFitConfig, NuFitResults);
void ProcessResults(std::vector<NuFitData*>, NuFitPDFs*, const NuFitConfig, std::vector<NuFitResults>);
void fitToFile(std::ofstream&, NuFitResults, const NuFitConfig);
void plotToFile(TFile*, NuFitData*, NuFitPDFs*, const NuFitConfig, NuFitResults);
std::vector<double> toCpdPerkton(std::vector<double>, const NuFitConfig, NuFitResults);

}  // namespace NuFitter


#endif
