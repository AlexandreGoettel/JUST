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
	std::vector<std::vector<double>> pcov_, std::vector<double> eff_,
	int errflg, int errflg_cov, const NuFitConfig config_, std::vector<std::vector<paramData>> vec,
	std::vector<std::vector<paramData>> vec_fixed) {

	popt = popt_;
	pcov = pcov_;
	efficiencies = eff_;
	errorflag = errflg;
	errorflag_cov = errflg_cov;
	paramVector = combineParamVectors(vec, vec_fixed, config_);
}

// @brief combine paramVectors to include free/constrained and fixed parameters
auto NuFitResults::combineParamVectors(std::vector<std::vector<paramData>> vec,
	std::vector<std::vector<paramData>> vec_fixed,
	const NuFitConfig config) -> std::vector<std::vector<paramData>> {
	// Insert the fixed params where they first appear in specieslist
	std::vector<std::vector<paramData>> output;
	auto iFree {0U}, iFixed {0U};
	while (true) {
		// Initialise to a high value to avoid infinite loops
		auto idx_pdf_free = config.npdfs;
		auto idx_pdf_fixed = config.npdfs;
		if (iFree < vec.size()) {
			idx_pdf_free = vec[iFree][0].idx_pdf;
		}
		if (iFixed < vec_fixed.size()) {
			idx_pdf_fixed = vec_fixed[iFixed][0].idx_pdf;
		}

		// Exit condition
		if (idx_pdf_free == config.npdfs && idx_pdf_fixed == config.npdfs)
			break;

		// Add the next pdf to the parameter list
		if (idx_pdf_free < idx_pdf_fixed) {
			output.push_back(vec[iFree]);
			iFree++;
		} else if (idx_pdf_free > idx_pdf_fixed){
			output.push_back(vec_fixed[iFixed]);
			iFixed++;
		} else {
			throw std::invalid_argument("A parameter should not be free and fixed simultaneously");
		}
	}

	return output;
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

}

}  // namespace NuFitter
