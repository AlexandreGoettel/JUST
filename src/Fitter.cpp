//! @file	include/Fitter.cpp
//! @brief	Implementation of class to fit MC PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard includes
#include <memory>
#include <vector>
#include <algorithm>
#include <cassert>  // Disabled if NDEBUG is defed
#include <iostream>
// ROOT includes
#include "Math/Util.h"
#include "TMath.h"
#include "TMatrixDSym.h"
// Project includes
#include "Fitter.h"


namespace NuFitter {
namespace MCFit {

NuFitContainer fitCtnr;

// TODO: if needed, move this to a more accessible location?
// @brief Multiply all vector elements by the same number
template <class Q, class P>
auto MultiplyVectorByScalar(std::vector<Q> &v, P k) -> void {
    transform(v.begin(), v.end(), v.begin(), [k](auto &c){ return c*k; });
}

// @brief Get the index of element el in vector v (first occurence)
template <class T>
auto getIndexOf(T el, std::vector<T> v) -> int {
	auto it = std::find(v.begin(), v.end(), el);

	// If element was found, calculate the index
	if (it != v.end()) {
		return it - v.begin();
	}
	else {
		return -1;
	}
}

// @brief Evaluate whether the bin between lower_edge and upper edge is
// inside of the bin range
template <class T>
auto NuFitContainer::InFitRange(T lower_edge, T upper_edge) -> bool {
    return upper_edge > config.emin && lower_edge < config.emax;
}

// @brief Constructor for NuFitContainer
// @brief Create new data/pdf vector objects with applied fit range cuts
// @brief Then separate free/fixed parameters and create the index maps
NuFitContainer::NuFitContainer(NuFitData *data_, NuFitPDFs *pdfs_,
	                           const NuFitConfig config_) {
	data = data_;
	pdfs = pdfs_;
	config = config_;

	// Adjust the data, pdf vectors according to config, save efficiencies
	// 1. Create the data vector for all histograms
    for (auto n = 0U; n < data->data.size(); n++) {
        auto bin_edges = data->bin_edges[n];
        auto hist_data = data->data[n];
        std::vector<double> current_hist;

        for (auto i = 0U; i < bin_edges.size()-1; i++) {
    		// Assuming bin_edges is ordered
    		// Fill data_vector with the hist data between emin and emax
    		if (InFitRange(bin_edges[i], bin_edges[i+1])) {
    			current_hist.push_back(hist_data[i]);
    		}
    	}
        data_vector.push_back(current_hist);
    }

	// 2. Create a vector of PDF vectors. Save efficiencies
    for (auto n = 0U; n < pdfs->pdfs.size(); n++) {
        std::vector<double> current_pdf;
        auto raw_pdf = pdfs->pdfs[n];
        auto current_efficiency {0.};
        auto bin_edges = pdfs->bin_edges[config.hist_id[n]-1];

        // Assuming bin_edges is orderd, fill pdfs between emin and emax
        for (auto i = 0U; i < bin_edges.size()-1; i++) {
			if (InFitRange(bin_edges[i], bin_edges[i+1])) {
				current_pdf.push_back(raw_pdf[i]);
			} else {
				current_efficiency += raw_pdf[i];
			}
		}

        // Renormalise vector
		MultiplyVectorByScalar(current_pdf, 1. / (1. - current_efficiency));

        // Save output
		pdf_vectors.push_back(current_pdf);
		efficiencies.push_back(1. - current_efficiency);
    }

    // Quick check that everything went according to plan
    for (auto i = 0U; i < config.hist_id.size(); i++) {
        auto n = config.hist_id[i];
        assert(data_vector[n-1].size() == pdf_vectors[i].size());
        assert(data_vector[n-1].size() != 0);
    }

	// 4. Fill the paramVector object which contains useful information
	//    about the parameters to be used in the fit
	std::vector<TString> used_names;
	for (auto i = 0U; i < fitCtnr.config.param_names.size(); i++) {
		paramData current_paramData {i, fitCtnr.config.hist_id[i]};

		auto name = fitCtnr.config.param_names[i];
		if (std::find(used_names.begin(), used_names.end(), name) == used_names.end()) {
			std::vector<paramData> tmp_paramVector;
			tmp_paramVector.push_back(current_paramData);
			paramVector.push_back(tmp_paramVector);
			used_names.push_back(name);
		} else {
			auto idx_name = getIndexOf(name, used_names);
			assert(idx_name != -1);
			paramVector[idx_name].push_back(current_paramData);
		}
	}
    n_params = paramVector.size();

    // 5. Create a vector of indices to map free params to pdfs
	// for (auto i = 0U; i < config.npdfs; i++) {
	// 	auto isFixed = config.param_fixed[i];
	// 	assert(isFixed == 0 || isFixed == 1 || isFixed == 2);
    //
	// 	if (isFixed != 1) {
	// 		idx_map.push_back(i);
	// 	} else {
	// 		idx_map_fixed.push_back(i);
	// 	}
	// }
	// n_fixed = idx_map_fixed.size();

    // 6. Make sure params are fixed/constr properly in the new setup
    // 7. If same var same pdf on same hist->raise error?

}

// Maybe return answer as a vector?
auto NuFitContainer::fitFunction(unsigned int npar, const double *par)
	                             -> std::vector<std::vector<double>> {
	// TODO: fixed params
	// Initialise output vector
	std::vector<std::vector<double>> fitFuncVal;
	for (auto el : data_vector) {
		std::vector<double> tmp (el.size(), 0);
		fitFuncVal.push_back(tmp);
	}

	// Fill the output vector with function values
	for (auto k = 0U; k < paramVector.size(); k++) {  // each parameter
		auto parDataVec = paramVector[k];
		for (auto h = 0U; h < parDataVec.size(); h++) {  // each hist in parData
			auto parData = parDataVec[h];
			for (auto j = 0U; j < data_vector[parData.idx_hist-1].size(); j++) {  // each bin
                fitFuncVal[parData.idx_hist-1][j] += pdf_vectors[parData.idx_pdf][j]
                    * par[k] * config.param_eff[parData.idx_pdf];
			}
		}
	}

	return fitFuncVal;
}
//
// @brief The fit function (parameters * pdfs)
// auto NuFitContainer::fitFunction(unsigned int i, unsigned int npar,
// 	                             const double *par, int nHist) -> double {
// 	// Add contributions from the free parameters
// 	auto yi {0.};
// 	for (auto j = 0U; j < npar; j++) {
// 		auto idx = fitCtnr.idx_map[j];  // Convert to param index space
// 		if (config.hist_id[idx] != nHist) continue;  // Ignore other hists
// 		yi += par[j] * pdf_vectors[idx][i];
// 	}
//
// 	// Now contributions from the fixed parameters
// 	// TODO: can speed-up with look-up since the result from below is constant
// 	for (auto j = 0U; j < fitCtnr.n_fixed; j++) {
// 		auto idx = fitCtnr.idx_map_fixed[j];  // Convert to param index space
// 		if (config.hist_id[idx] != nHist) continue;  // Ignore other hists
// 		yi += fitCtnr.config.param_initial_guess[idx] * pdf_vectors[idx][i];
// 	}
//
// 	return yi;
// }

// @brief Define Minuit-Style binned poisson likelihood (extended)
auto NuFitContainer::NLL_poisson(int npar, const double *par) -> double {
	// Following Baker&Cousins 1983 definition on page 439
	auto chi_sqr_lambda_p { 0. };
	auto fitFuncVal = fitFunction(npar, par);

	// Make sure data and model shapes are compatible
	assert(data_vector.size() == fitFuncVal.size());
	assert(data_vector[0].size() == fitFuncVal[0].size());

	// Loop over data
	for (auto i = 0U; i < data_vector.size(); i++) {
		for (auto j = 0U; j < data_vector[i].size(); j++) {
			auto yi = fitFuncVal[i][j];
			auto ni = data_vector[i][j];

			chi_sqr_lambda_p += yi - ni + ni*(ROOT::Math::Util::EvalLog(ni) -
			                                  ROOT::Math::Util::EvalLog(yi));
		}
	}

	return chi_sqr_lambda_p;
}

// @brief Define MUST Likelihood (for comparison purposes)
// TODO: Re-introduce
// auto NuFitContainer::NLL_MUST(int npar, const double *par) -> double {
// 	auto logL { 0. };
//
// 	// Loop over possible histograms
// 	for (auto n = 0U; n < data_vector.size(); n++) {
// 		auto data_hist = data_vector[n];
// 		auto nbins = data_hist.size();
// 		for (auto i = 0U; i < nbins; i++) {
// 			auto yi = fitFunction(i, npar, par, n);
// 			auto ni = data_hist[i];
//
// 			if (ni < 0)
// 				continue;
// 			else if (ni == 0)
// 				logL += -yi;
// 			else
// 				logL += ni*ROOT::Math::Util::EvalLog(yi)-yi-TMath::LnGamma(ni+1);
// 		}
// 	}
// 	return -logL;
// }

// @brief MinuitManager constructor
MinuitManager::MinuitManager(const NuFitConfig config_) {
	config = config_;
	errorflag = 0;
}

// @brief Initialise Minuit with params etc.
auto MinuitManager::initMinuit() -> void {
	// Create a new instance
	gMinuit = new TMinuit();  // TODO: Add protection -> v0.3?
	gMinuit->SetFCN(fcn);
	// Set each parameter in Minuit
	for (auto j = 0U; j < fitCtnr.n_params; j++) {
		// Convert index space
		auto i = fitCtnr.paramVector[j][0].idx_pdf;
		// Give the parameter information to Minuit
		gMinuit->mnparm(j, config.param_names[i],
			config.param_initial_guess[i], config.param_stepsize[i],
			config.param_lowerlim[i], config.param_upperlim[i], errorflag);
	}
}

// @brief start the minimization process by executing Minuit commands
auto MinuitManager::callMinuit() -> void {
	// Give Minuit a list of commands through arglist
	double arglist[2];

	// We are doing maximum likelihood fits: errors at +0.5 logL
	arglist[0] = 0.5;
	gMinuit ->mnexcm("SET ERR", arglist, 1, errorflag);

	// STR=1: standard fit
	// STR=2: additional search around found minimum, needs derivatives
	arglist[0] = 2;
	gMinuit->mnexcm("SET STR", arglist, 1, errorflag);

	// Call MIGRAD (+ SIMPLEX method if Migrad fails)
	arglist[0] = 50000;
	arglist[1] = 0.001;
	gMinuit->mnexcm("MINIMIZE", arglist, 2, errorflag);

	// Optional: call extra Hesse calculation
	if (config.doHesse) {
		gMinuit->mnexcm("HESSE", arglist, 2, errorflag);
	}
	// Optional: do exact non-linear error calculation
	if (config.doMinos) {
		gMinuit->mnexcm("MINOS", arglist, 2, errorflag);
	}
}

// @brief Convert fit results to vectors, store in member variables popt/pcov
auto MinuitManager::getResults() -> NuFitResults {
	// Initialise variables
	double x, _;
    unsigned int nparams = fitCtnr.n_params;
	std::vector<double> popt(nparams);
	std::vector<std::vector<double>> pcov(nparams,
	                                      std::vector<double>(nparams));

	// Fill parameter vector
    for (auto i = 0U; i < fitCtnr.n_params; i++) {
        auto j = fitCtnr.paramVector[i][0].idx_pdf;
		gMinuit->GetParameter(i, x, _);
		popt[j] = x;
	}
	// Fill fixed parameter value where estimate would be
	// TODO

	// Get the covariance matrix
	double mat[nparams][nparams];
	gMinuit->mnemat(&mat[0][0], nparams);
	std::vector<std::vector<double>> pcov_(nparams,
	                                       std::vector<double>(nparams));
	// Convert to vector<vector>
	for (auto i = 0U; i < nparams; i++) {
		for (auto j = 0U; j < nparams; j++) {
			pcov_[i][j] = mat[i][j];
		}
	}

	// Expand pcov to include fixed params as well
    // TODO
	// For each row, this var gives the idx of free params (inverse of idx_map)
	// auto idx_inverse {0U};
	// for (auto iRow = 0U; iRow < config.npdfs; iRow++){
	// 	// Construct row
	// 	std::vector<double> this_row(config.npdfs);
	// 	// If iRow is fixed, row[iRow]=1
	// 	if (std::find(fitCtnr.idx_map_fixed.begin(), fitCtnr.idx_map_fixed.end(),
	// 	              iRow) == fitCtnr.idx_map_fixed.end()) {
	// 		// If not fixed, get data from pcov_
	// 		for (auto i = 0U; i < fitCtnr.idx_map.size(); i++) {
	// 			auto j = fitCtnr.idx_map[i];
	// 			this_row[j] = pcov_[idx_inverse][i];
	// 		}
	// 		idx_inverse += 1;
	// 	}
	// 	pcov[iRow] = this_row;
	// }

    // Get the status of the covariance matrix calculation
    int tmp_;
    gMinuit->mnstat(_, _, _, tmp_, tmp_, errorflag_cov);

	auto results = NuFitResults(popt, pcov, fitCtnr.efficiencies,
                                errorflag, errorflag_cov);
	return results;
}

// @brief Used by TMinuit to sample the likelihood
auto fcn(int &npar, double *gin, double &f, double *par, int iflag) -> void {
	// Calculate the log-likelihood according to different methods
	if (fitCtnr.config.likelihood.compare("poisson") == 0) {
		f = fitCtnr.NLL_poisson(npar, par);
	// } else if (fitCtnr.config.likelihood.compare("must") == 0) {
	// 	f = fitCtnr.NLL_MUST(npar, par);
	} else {
		throw std::invalid_argument("'" + fitCtnr.config.likelihood +
			"' is not a valid likelihood.\nAllowed: ['poisson']" +
			"\nPay attention to the capitalisation!");
	}

	// In case of parameter constraints, add (Gaussian) pull terms here
	// TODO
	// for (auto i = 0U; i < fitCtnr.n_params; i++) {
	// 	auto j = fitCtnr.idx_map[i];  // Convert to free param index space
	// 	if (fitCtnr.config.param_fixed[j] != 2) continue;
	//
	// 	auto diff = par[i] - fitCtnr.config.param_initial_guess[j];
	// 	auto sigma = fitCtnr.config.param_stepsize[j];
	// 	f += 0.5*diff*diff/sigma/sigma;
	// }
}

// @brief Perform a binned likelihood fit of the pdfs on the data
// @param data data to be fitted
// @param pdfs MC PDFs to fit to
// @param config container for fit options / variables
// @return NuFitResults object containing relevant fit results info
auto Fit(NuFitData *data, NuFitPDFs *pdfs, const NuFitConfig config)
		-> NuFitResults {
	// 1. Create the NuFitContainer object -> formats data according to config
	// 1.1 Overwrite NuFitContainer at static location in MCFit scope
	// TODO: Need memory leak protection?
	new (&fitCtnr) NuFitContainer(data, pdfs, config);

	// 2. Prepare TMinuit
	auto manager = std::make_unique<MinuitManager>(config);
	manager->initMinuit();

	// 3. Start minimization
	manager->callMinuit();

	// 4. Return results
	return manager->getResults();
}

// @brief Fit for each entry in vector<NuFitData*>
// @param data vector of NuFitData pointers
// @param pdfs pointer to NuFitPDFs with the MC PDFs
// @param config pointer to the fit config variables
// @return vector of NuFitResults*, one for each NuFitData* in data
auto Fit(std::vector<NuFitData*> data, NuFitPDFs *pdfs,
	     const NuFitConfig config) -> NuFitResults {
	auto results = Fit(data[0], pdfs, config);
	for (auto i = 1u; i < data.size(); i++) {
		results.addResults(Fit(data[i], pdfs, config));
	}
	return results;
}

}  // namespace MCFit
}  // namespace NuFitter
