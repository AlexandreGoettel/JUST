//! @file	src/FitResults.cpp
//! @brief	Implementation of class to contain fit results
//! @author	Alexandre Göttel
//! @date	2021-04-13

// includes
#include "FitResults.h"

namespace NuFitter {

// @brief Constructor
NuFitResults::NuFitResults(std::vector<double> popt_,
	                       std::vector<std::vector<double>> pcov_) {
   popt = popt_;
   pcov = pcov_;
}

// @brief Add results to existing NuFitResults. Used for ToyData fits.
auto NuFitResults::addResults(NuFitResults results) -> void {
	// TODO
}

}  // namespace NuFitter
