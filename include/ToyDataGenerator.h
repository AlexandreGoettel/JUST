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
	NuFitToyData(const NuFitConfig, NuFitPDFs*);  // constructor
	~NuFitToyData() = default;  // destructor
	NuFitToyData(const NuFitToyData&) = default;  // copy constructor
	NuFitToyData(NuFitToyData&&) = default;  // move constructor
	NuFitToyData &operator=(const NuFitToyData&) = default;  // copy assignment
	NuFitToyData &operator=(NuFitToyData&&) = default;  // move assignment

public: // Variables
	NuFitPDFs *pdfs;
	void loadDataset(unsigned int);
	NuFitData *dataset = nullptr;

private:
	NuFitConfig config;
	std::vector<TH1D*> histogr;
	std::vector<std::vector<double>> bin_edges;
	std::vector<unsigned int> hist_ids;
};

namespace ToyData {
	NuFitToyData *Initialise(const NuFitConfig, NuFitPDFs*);
}  // namespace ToyData
}  // namespace NuFitter


#endif
