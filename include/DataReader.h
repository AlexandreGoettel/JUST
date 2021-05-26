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

//============================================================================
// Class definition

namespace NuFitter {

// Forward declaration
class NuFitConfig;

class NuFitData {
public:  // Constructros and operator assigments
	NuFitData(std::vector<std::vector<double>>, std::vector<std::vector<double>>,
	          std::vector<TH1D*>);  // constructor
	~NuFitData() = default;  // destructor
	NuFitData(const NuFitData&) = default;  // copy constructor
	NuFitData(NuFitData&&) = default;  // move constructor
	NuFitData &operator=(const NuFitData&) = default;  // copy assignment
	NuFitData &operator=(NuFitData&&) = default;  // move assignment

public: // Variables
	std::vector<std::vector<double>> data, bin_edges;
	std::vector<TH1D*> data_histograms;
};

class NuFitPDFs {
public:  // Constructors and operator assigments
	NuFitPDFs(std::vector<std::vector<double>>, std::vector<double>,
	          std::vector<TH1D*>);  // constructor
	~NuFitPDFs() = default;  // destructor
	NuFitPDFs(const NuFitPDFs&) = default;  // copy constructor
	NuFitPDFs(NuFitPDFs&&) = default;  // move constructor
	NuFitPDFs &operator=(const NuFitPDFs&) = default;  // copy assignment
	NuFitPDFs &operator=(NuFitPDFs&&) = default;  // move assignment

public:  // Variables
	std::vector<double> bin_edges;
	std::vector<std::vector<double>> pdfs;
	std::vector<TH1D*> pdf_histograms;
};

namespace Data {
	NuFitData *Read(const NuFitConfig config);
}  // namespace CMDLParser

namespace PDFs {
	NuFitPDFs* Read(const NuFitConfig config);
}  // namespace ConfigParser

}  // namespace NuFitter


#endif
