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
#include "ToyDataGenerator.h"

//============================================================================
// Method definitions

namespace NuFitter {

// TODO: maybe a class here would be useful after all?
// TODO: a lot of these refs can/should be const refs
void ProcessResults(NuFitData*&, NuFitPDFs*&, NuFitResults*&);
void ProcessResults(NuFitToyData*&,NuFitPDFs*&, NuFitPDFs*&, std::vector<NuFitResults*>&);

void fitToFile(std::ofstream&, NuFitResults*&);
void plotToFile(TFile*&, NuFitData*&, NuFitPDFs*&, NuFitResults*&);
void DrawPDFs(TFile*&, NuFitPDFs*&, NuFitPDFs*&, NuFitResults*&);
void createParamTree(TFile*&, NuFitResults*&);
std::vector<double> toCpdPerkton(std::vector<double>, NuFitResults*&);

}  // namespace NuFitter


#endif
