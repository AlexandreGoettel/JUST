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

// @brief Copy cov matrix to corr matrix, output
auto NuFitResults::getCorrMatrix() -> std::vector<std::vector<double>> {
	auto n = pcov.size();
	std::vector<std::vector<double>> out(n, std::vector<double>(n));
	for (auto i = 0U; i < n; i++) {
		auto var_i = pcov[i][i];
		for (auto j = 0U; j < n; j++) {
			auto var_j = pcov[j][j];
			if (var_i == 0 && var_j == 0) {
				out[i][j] = 1;
			} else if (var_i*var_j != 0){
				out[i][j] = pcov[i][j] / std::sqrt(var_i*var_j);
			}
		}
	}
	return out;
}

// @brief Add results to existing NuFitResults. Used for ToyData fits.
auto NuFitResults::addResults(NuFitResults results) -> void {
	// TODO
}

}  // namespace NuFitter
