//! @file	src/FitResults.cpp
//! @brief	Implementation of class to contain fit results
//! @author	Alexandre Göttel
//! @date	2021-04-13

// includes
#include <cmath>
#include <stdexcept>
#include "FitResults.h"

namespace NuFitter {

// @brief Constructor
NuFitResults::NuFitResults(const NuFitConfig config_, std::vector<double> popt_,
	std::vector<std::vector<double>> pcov_, std::vector<double> eff_,
	int errflg, int errflg_cov,
	std::vector<std::vector<paramData>> vec,
    std::vector<std::vector<paramData>> vec_fixed) {
	// Initialise variables
	config = config_;
	popt = popt_;
	pcov = pcov_;
	efficiencies = eff_;
	errorflag = errflg;
	errorflag_cov = errflg_cov;

	// Combine paramVectors for the OutputManager
	combineParamVectors(vec, vec_fixed);
}

// @brief combine free and fixed parameter info in one paramVector
auto NuFitResults::combineParamVectors(std::vector<std::vector<paramData>> vec,
	std::vector<std::vector<paramData>> vec_fixed) -> void {

	std::vector<double> popt_extended;

	// Fill paramVector in the order that the parameters were seen by parser
	auto iFree {0U}, iFixed {0U};
	for (auto i = 0U; i < vec.size()+vec_fixed.size(); i++) {
		auto idx_pdf_free = config.npdfs;  // Initialise to a larger number
		auto idx_pdf_fixed = config.npdfs;  // For the comparison below
		if (iFree < vec.size()) {
			idx_pdf_free = vec[iFree][0].idx_pdf;
		}
		if (iFixed < vec_fixed.size()) {
			idx_pdf_fixed = vec_fixed[iFixed][0].idx_pdf;
		}

		if (idx_pdf_free < idx_pdf_fixed) {
			paramVector.push_back(vec[iFree]);
			popt_extended.push_back(popt[iFree]);
			iFree++;
		} else if (idx_pdf_free > idx_pdf_fixed){
			paramVector.push_back(vec_fixed[iFixed]);
			popt_extended.push_back(config.param_initial_guess[idx_pdf_fixed]);
			iFixed++;
		} else {
			throw std::range_error("One variable was seen as free and fixed simultaneously");
		}
	}

	// Update popt to include free parameters
	popt = popt_extended;
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
