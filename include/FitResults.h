//! @file	include/NuFitResults.h
//! @brief	Declaration of class to contain fit results
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef NuFitResults_H_
#define NuFitResults_H_

//============================================================================
// Includes
#include <vector>

//============================================================================
// Class definition

namespace NuFitter {

class NuFitResults {
// Constructors and operator assignments
public:
	NuFitResults(std::vector<double>, std::vector<std::vector<double>>,
		         std::vector<double>);  // constructor
	~NuFitResults() = default;  // destructor
	NuFitResults(const NuFitResults&) = default;  // copy constructor
	NuFitResults(NuFitResults&&) = default;  // move constructor
	NuFitResults &operator=(const NuFitResults&) = default;  // copy assignment
	NuFitResults &operator=(NuFitResults&&) = default;  // move assignment

public: // Member variables
	std::vector<double> popt, efficiencies;
	std::vector<std::vector<double>> pcov;

public: // Functions
	void addResults(NuFitResults);
	std::vector<double> getUncertainties();
	std::vector<std::vector<double>> getCorrMatrix();
};

}  // namespace NuFitter


#endif
