//! @file	include/NuFitResults.h
//! @brief	Declaration of class to contain fit results
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef NuFitResults_H_
#define NuFitResults_H_

//============================================================================
// Includes
#include <vector>
#include "Parser.h"

//============================================================================
// Class definition

namespace NuFitter {

class NuFitResults {
public: // Constructors and operator assignments
	NuFitResults(std::vector<double>, std::vector<std::vector<double>>,
		         std::vector<double>, int, int,
	             std::vector<std::vector<paramData>>,
			 	 std::vector<std::vector<paramData>>, double);  // constructor
	~NuFitResults();  // destructor
	NuFitResults(const NuFitResults&) = delete;  // copy constructor
	NuFitResults(NuFitResults&&) = delete;  // move constructor
	NuFitResults &operator=(const NuFitResults&) = delete;  // copy assignment
	NuFitResults &operator=(NuFitResults&&) = delete;  // move assignment

public: // Member variables
	std::vector<double> popt, efficiencies;
	std::vector<std::vector<double>> pcov;
	unsigned int errorflag, errorflag_cov;
	double chi_sqr_ndof;
	std::vector<std::vector<paramData>> paramVector;

public: // Functions
	std::vector<double> getUncertainties();
	std::vector<std::vector<double>> getCorrMatrix();

private:
	std::vector<std::vector<paramData>> combineParamVectors(
		std::vector<std::vector<paramData>>,
		std::vector<std::vector<paramData>>);
};
}  // namespace NuFitter

#endif
