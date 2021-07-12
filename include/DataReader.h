//! @file	include/DataReader.h
//! @brief	Declaration of classes to read data and PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef DATAREADER_H_
#define DATAREADER_H_

//============================================================================
// Includes
#include <vector>
#include "TH1.h"
#include "TFile.h"

//============================================================================
// Class definition

namespace NuFitter {

class NuFitData {
public:  // Constructors and operator assigments
	NuFitData(std::vector<std::vector<double>>, std::vector<std::vector<double>>,
	          std::vector<TH1D*>, std::vector<unsigned int>);  // constructor
	NuFitData(std::vector<std::vector<double>>, std::vector<std::vector<double>>,
	          std::vector<unsigned int>);  // constructor
	NuFitData();
	~NuFitData();  // destructor
	NuFitData(const NuFitData&) = delete;  // copy constructor
	NuFitData(NuFitData&&) = delete;  // move constructor
	NuFitData &operator=(const NuFitData&) = delete;  // copy assignment
	NuFitData &operator=(NuFitData&&) = delete;  // move assignment

public: // Variables
	std::vector<std::vector<double>> data, bin_edges;
	std::vector<TH1D*> data_histograms;
	std::vector<unsigned int> hist_ids;
	void Read();
private:
	TFile *file_data = nullptr;
};

class NuFitPDFs {
public:  // Constructors and operator assigments
	NuFitPDFs();  // constructor
	~NuFitPDFs();  // destructor
	NuFitPDFs(const NuFitPDFs&) = delete;  // copy constructor
	NuFitPDFs(NuFitPDFs&&) = delete;  // move constructor
	NuFitPDFs &operator=(const NuFitPDFs&) = delete;  // copy assignment
	NuFitPDFs &operator=(NuFitPDFs&&) = delete;  // move assignment

public:  // Variables
	std::vector<std::vector<double>> pdfs, bin_edges;
	std::vector<TH1D*> pdf_histograms;
	void Read(const std::vector<TString>&, const unsigned int&, const std::vector<unsigned int>&);
private:
	TFile *file_pdf = nullptr;
};
}  // namespace NuFitter


#endif
