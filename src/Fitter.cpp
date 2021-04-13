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

// @brief Define Minuit-Style binned poisson likelihood
auto NuFitContainer::NLL(int npar, const double *par) -> double {
	// Following Baker&Cousins 1983 definition on page 439
	auto nbins = data_vector.size();
	auto chi_sqr_lambda_p { 0. };
	for (auto i = 0U; i < nbins; i++) {
		auto yi = fitFunction(i, npar, par);
		// auto yi = par[0]*signal_vector[i] + par[1]*background_vector[i];
		auto ni = data_vector[i];

		chi_sqr_lambda_p += yi - ni + ni*(ROOT::Math::Util::EvalLog(ni) -
		                                  ROOT::Math::Util::EvalLog(yi));
	}
	return chi_sqr_lambda_p;
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
		gMinuit->mnparm(i, config.param_names[i],
			config.param_initial_guess[i], config.param_stepsize[i],
			config.param_lowerlim[i], config.param_upperlim[i], errorflag);
	}
}

// @brief start the minimization process by executing Minuit commands
auto MinuitManager::callMinuit() -> void {
	double arglist[2];
	arglist[0] = 1e-15L;
	gMinuit->mnexcm("SET EPS", arglist, 1, errorflag);

	arglist[0] = 0.5;
	gMinuit ->mnexcm("SET ERR", arglist, 1, errorflag);

	arglist[0] = 2;
	gMinuit->mnexcm("SET STR", arglist, 1, errorflag);

	arglist[0] = 50000;
	arglist[1] = 0.001;
	gMinuit->mnexcm("CALL FCN", arglist, 2, errorflag);
	gMinuit->mnexcm("MINIMIZE", arglist, 2, errorflag);

	// TODO: Use config to decide whether to call hesse/minos?
	if (config.doHesse)
		gMinuit->mnexcm("HESSE", arglist, 2, errorflag);
	if (config.doMinos)
		gMinuit->mnexcm("MINOS", arglist, 2, errorflag);
}

// @brief Perform a binned likelihood fit of the pdfs on the data
// @param data data to be fitted
// @param pdfs MC PDFs to fit to
// @param config container for fit options / variables
// @return NuFitResults object containing relevant fit results info
auto Fit(NuFitData *data, NuFitPDFs *pdfs, const NuFitConfig config)
		-> NuFitResults* {
	// 1. Create the NuFitContainer object -> formats data according to config
	// auto fitCtnr = std::make_unique<NuFitContainer>(data, pdfs, config);
	// 1. Overwrite NuFitContainer at static location
	new (&fitCtnr) NuFitContainer(data, pdfs, config);

	// 2. Prepare TMinuit
	auto manager = std::make_unique<MinuitManager>(config);
	manager->initMinuit();

	// 3. Start minimization
	manager->callMinuit();

	// 4. Pass output to NuFitResults
	// TODO
	auto results = std::make_unique<NuFitResults>();
	return results.get();
}

// @brief Fit for each entry in vector<NuFitData*>
// @param data vector of NuFitData pointers
// @param pdfs pointer to NuFitPDFs with the MC PDFs
// @param config pointer to the fit config variables
// @return vector of NuFitResults*, one for each NuFitData* in data
auto Fit(std::vector<NuFitData*> data, NuFitPDFs *pdfs,
	     const NuFitConfig config) -> NuFitResults* {
	auto results = Fit(data[0], pdfs, config);
	for (auto i = 1u; i < data.size(); i++) {
		results->addResults(Fit(data[i], pdfs, config));
	}
	return results;
}

// @brief Used by TMinuit to sample the likelihood
auto fcn(int &npar, double *gin, double &f, double *par, int iflag) -> void {
	f = fitCtnr.NLL(npar, par);
}

}  // namespace MCFit
}  // namespace NuFitter
