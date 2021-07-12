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
#include <cstring>
// ROOT includes
#include "TFile.h"
// Project includes
#include "FitResults.h"
#include "DataReader.h"
#include "ToyDataGenerator.h"

//============================================================================
// Method definitions

namespace NuFitter {

class OutputManager {
public:  // Constructors and operator assigments
	// OutputManager(NuFitData*&, NuFitPDFs*&, NuFitResults*&);
	// OutputManager(NuFitToyData*&,NuFitPDFs*&, NuFitPDFs*&, std::vector<NuFitResults*>&);
	OutputManager() = default;
	~OutputManager();  // destructor
	OutputManager(const OutputManager&) = default;  // copy constructor
	OutputManager(OutputManager&&) = default;  // move constructor
	OutputManager &operator=(const OutputManager&) = default;  // copy assignment
	OutputManager &operator=(OutputManager&&) = default;  // move assignment

public:
	void initRootFile(std::string&);
	void closeRootFile();
	void writeParamTree(NuFitResults*&);
	void writeFitTree(std::vector<NuFitResults*>&);
	void writeToyTree();
	void plotToFile(NuFitData*&, NuFitPDFs*&, NuFitResults*&);
	void drawPDFs(NuFitPDFs*&, NuFitPDFs*&);
	void fitToFile(std::ofstream&, NuFitResults*&);

private:
	void makePDFsSum(NuFitData*&, NuFitPDFs*&, NuFitResults*&);

private:
	TFile *f;
	std::vector<TH1D*> *PDFsSum = nullptr;
	std::vector<std::vector<paramData>> paramVector;
};

// TODO: a lot of these refs can/should be const refs
void ProcessResults(NuFitData*&, NuFitPDFs*&, NuFitResults*&);
void ProcessResults(NuFitToyData*&,NuFitPDFs*&, NuFitPDFs*&, std::vector<NuFitResults*>&);

std::vector<double> toCpdPerkton(std::vector<double>, NuFitResults*&);

}  // namespace NuFitter


#endif
