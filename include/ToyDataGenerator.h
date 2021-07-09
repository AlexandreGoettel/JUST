//! @file	include/ToyDataGenerator.h
//! @brief	Declaration of class to generate toy data
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef TOYDATAGENERATOR_H_
#define TOYDATAGENERATOR_H_

//============================================================================
// Standard includes
#include <vector>
// ROOT includes
#include "TH1D.h"
// Project includes
#include "Parser.h"
#include "DataReader.h"

//============================================================================
// Class definition

namespace NuFitter {

// Forward declaration
class NuFitData;

class NuFitToyData {
public:  // Constructros and operator assigments
	NuFitToyData(NuFitPDFs*&);  // constructor
	~NuFitToyData();  // destructor
	NuFitToyData(const NuFitToyData&) = delete;  // copy constructor
	NuFitToyData(NuFitToyData&&) = delete;  // move constructor
	NuFitToyData &operator=(const NuFitToyData&) = delete;  // copy assignment
	NuFitToyData &operator=(NuFitToyData&&) = delete;  // move assignment

public:
	void loadDataset(unsigned int);
	NuFitPDFs *pdfs;
	NuFitData *dataset = nullptr;
	TH1D *hdata = nullptr;

private:
	std::vector<TH1D*> histogr;
	std::vector<std::vector<double>> bin_edges;
	std::vector<unsigned int> hist_ids;
};

namespace ToyData {
	NuFitToyData *Initialise(NuFitPDFs*&);
}  // namespace ToyData
}  // namespace NuFitter


#endif
