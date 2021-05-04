//! @file	include/Fitter.cpp
//! @brief	Implementation of class to fit MC PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard includes
#include <memory>
#include <vector>
#include <cassert>  // Disabled if NDEBUG is defed
#include <iostream>
// ROOT includes
#include "Math/Util.h"
#include "TMatrixDSym.h"
// Project includes
#include "Fitter.h"


namespace NuFitter {
namespace MCFit {

NuFitContainer fitCtnr;

// @brief Constructor for NuFitContainer
NuFitContainer::NuFitContainer(NuFitData *data_, NuFitPDFs *pdfs_,
	                           const NuFitConfig config_) {
	data = data_;
	pdfs = pdfs_;
	config = config_;

	// Create the data, pdf vectors according to config, save efficiencies
	// 1. Get relevant variables
	auto data_raw = data->data;
	auto bin_edges = data->bin_edges;
	assert(bin_edges == pdfs->bin_edges);

	// TODO: Extend to multiple histograms..
	auto emin = config.emin;
	auto emax = config.emax;

	// 2. Create the data vector
	for (auto i = 0U; i < bin_edges.size()-1; i++) {
		// Assuming the bin_edges is ordered
		// Fill data_vector with the raw_data between emin and emax
		if (bin_edges[i+1] > emin && bin_edges[i] < emax) {
			data_vector.push_back(data_raw[i]);
		}
	}

	// 3. Create a vector of PDF vectors. Save efficiencies
	for (auto el : pdfs->pdfs) {
		auto current_efficiency {0.};
		std::vector<double> current_pdf;

		for (auto i = 0U; i < bin_edges.size()-1; i++) {
			if (bin_edges[i+1] > emin && bin_edges[i] < emax) {
				current_pdf.push_back(el[i]);
			} else {
				current_efficiency += el[i];
			}
		}
		pdf_vectors.push_back(current_pdf);
		efficiencies.push_back(current_efficiency);
	}
	assert(data_vector.size() == pdf_vectors[0].size());
}

// @brief The fit function (parameters * pdfs)
auto NuFitContainer::fitFunction(unsigned int i, unsigned int npar, const double *par)
		-> double {
	auto yi {0.};
	for (auto j = 0U; j < npar; j++) {
		yi += par[j] * pdf_vectors[j][i];
	}
	return yi;
}

// @brief Define Minuit-Style binned poisson likelihood (extended)
auto NuFitContainer::NLL_extended(int npar, const double *par) -> double {
	// Following Baker&Cousins 1983 definition on page 439
	auto nbins = data_vector.size();
	auto chi_sqr_lambda_p { 0. };
	for (auto i = 0U; i < nbins; i++) {
		auto yi = fitFunction(i, npar, par);
		auto ni = data_vector[i];

		chi_sqr_lambda_p += yi - ni + ni*(ROOT::Math::Util::EvalLog(ni) -
		                                  ROOT::Math::Util::EvalLog(yi));
	}
	return chi_sqr_lambda_p;
}

// @brief Define standard binned poisson likelihood
auto NuFitContainer::NLL_poisson(int npar, const double *par) -> double {
	auto nbins = data_vector.size();
	auto logL { 0. };
	for (auto i = 0U; i < nbins; i++) {
		auto yi = fitFunction(i, npar, par);
		auto ni = data_vector[i];

		logL += ni*ROOT::Math::Util::EvalLog(yi) - yi;
	}
	return -logL;
}

// @brief MinuitManager constructor
MinuitManager::MinuitManager(const NuFitConfig config_) {
	config = config_;
	errorflag = 0;
}

// @brief Initialise Minuit with params etc.
auto MinuitManager::initMinuit() -> void {
	// Create a new instance
	gMinuit = new TMinuit();  // TODO: Add protection?
	gMinuit->SetFCN(fcn);
	// Set each parameter in Minuit
	for (auto i = 0U; i < config.nparams; i++) {
		// Give the parameter information to Minuit
		gMinuit->mnparm(i, config.param_names[i],
			config.param_initial_guess[i], config.param_stepsize[i],
			config.param_lowerlim[i], config.param_upperlim[i], errorflag);

		// Fix parameters
		if (config.param_fixed[i] == 1) {
			gMinuit->FixParameter(i);
		}
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
	gMinuit->mnexcm("CALL FCN", arglist, 2, errorflag);
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
auto MinuitManager::getResults() -> void {
	// Get the covariance matrix
	TMatrixDSym mat(config.nparams);
	gMinuit->mnemat(mat.GetMatrixArray(), config.nparams);

	// Convert to vector
	for (auto i = 0U; i < config.nparams; i++) {
		std::vector<double> pcov_row;
		for (auto j = 0U; j < config.nparams; j++) {
			pcov_row.push_back(mat[i][j]);
		}
		pcov.push_back(pcov_row);
	}

	// Put the parameter estimates into another vector
	double x, sigmax;
	for (auto i = 0U; i < config.nparams; i++) {
		gMinuit->GetParameter(i, x, sigmax);
		popt.push_back(x);
		// Could also get an uncertainty vector here.
		popt_err.push_back(sigmax);
	}
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
	manager->getResults();

	// 4. Pass output to NuFitResults
	// TODO
	auto results = NuFitResults(manager->popt, manager->popt_err, manager->pcov);
	return results;
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

// @brief Used by TMinuit to sample the likelihood
auto fcn(int &npar, double *gin, double &f, double *par, int iflag) -> void {
	if (fitCtnr.config.likelihood.compare("poisson") == 0) {
		f = fitCtnr.NLL_poisson(npar, par);
	} else if (fitCtnr.config.likelihood.compare("extended") == 0) {
		f = fitCtnr.NLL_extended(npar, par);
	} else {
		throw std::invalid_argument("'" + fitCtnr.config.likelihood +
			"' is not a valid likelihood.\nAllowed: ['poisson', 'extended']" +
			"\nPay attention to the capitalisation!");
	}
}

}  // namespace MCFit
}  // namespace NuFitter
