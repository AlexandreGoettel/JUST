//! @file	src/FitResults.cpp
//! @brief	Implementation of class to contain fit results
//! @author	Alexandre Göttel
//! @date	2021-04-13

// includes
#include <cmath>
#include "FitResults.h"

namespace NuFitter {

// @brief Constructor
NuFitResults::NuFitResults(std::vector<double> popt_,
	                       std::vector<std::vector<double>> pcov_,
						   std::vector<double> eff_) {
	popt = popt_;
	pcov = pcov_;
	efficiencies = eff_;
}

// @brief Extract the uncertainties from the correlation matrix
auto NuFitResults::getUncertainties() -> std::vector<double> {
	std::vector<double> popt_err;
	for (auto i = 0U; i < pcov.size(); i++) {
		popt_err.push_back(std::sqrt(pcov[i][i]));
	}
	return popt_err;
}

// @brief Add results to existing NuFitResults. Used for ToyData fits.
auto NuFitResults::addResults(NuFitResults results) -> void {
	// TODO
}

}  // namespace NuFitter
