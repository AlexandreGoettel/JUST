//! @file	include/NuFitResults.h
//! @brief	Declaration of class to contain fit results
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef NuFitResults_H_
#define NuFitResults_H_

//============================================================================
// Includes
#include <vector>
#include "TString.h"
#include "Parser.h"

//============================================================================
// Class definition

namespace NuFitter {

struct paramData {
	unsigned int idx_pdf;
	unsigned int idx_hist;
};

class NuFitResults {
public: // Constructors and operator assignments
	NuFitResults(const NuFitConfig, std::vector<double>, std::vector<std::vector<double>>,
	             std::vector<double>, int, int,
	             std::vector<std::vector<paramData>>,
	             std::vector<std::vector<paramData>>);  // constructor
	~NuFitResults() = default;  // destructor
	NuFitResults(const NuFitResults&) = default;  // copy constructor
	NuFitResults(NuFitResults&&) = default;  // move constructor
	NuFitResults &operator=(const NuFitResults&) = default;  // copy assignment
	NuFitResults &operator=(NuFitResults&&) = default;  // move assignment

public: // Member variables
	NuFitConfig config;
	std::vector<double> popt, efficiencies;
	std::vector<std::vector<double>> pcov;
	unsigned int errorflag, errorflag_cov;
	std::vector<std::vector<paramData>> paramVector;

public: // Functions
	void addResults(NuFitResults);
	std::vector<double> getUncertainties();
	std::vector<std::vector<double>> getCorrMatrix();

private:
	void combineParamVectors(std::vector<std::vector<paramData>>,
		                     std::vector<std::vector<paramData>>);
};

}  // namespace NuFitter


#endif
