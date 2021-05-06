//! @file	src/FitResults.cpp
//! @brief	Implementation of class to contain fit results
//! @author	Alexandre Göttel
//! @date	2021-04-13

// includes
#include "FitResults.h"

namespace NuFitter {

// @brief Constructor
NuFitResults::NuFitResults(std::vector<double> popt_,
	std::vector<double> popt_err_, std::vector<double> eff_) {
   popt = popt_;
   popt_err = popt_err_;
   efficiencies = eff_;
}

// @brief Add results to existing NuFitResults. Used for ToyData fits.
auto NuFitResults::addResults(NuFitResults results) -> void {
	// TODO
}

}  // namespace NuFitter
