//! @file	src/OutputManager.cpp
//! @brief	Implementation of class to process fit results
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Includes
#include <iostream>
#include "TFile.h"
#include "OutputManager.h"

namespace NuFitter {

// @brief For now, simply plot the results (simple fit example)
auto ProcessResults(NuFitData *data, NuFitPDFs *pdfs, const NuFitConfig config, NuFitResults results) -> void {
	// Open file to save the plots in
	TFile *f = new TFile(config.output_name.c_str(), "RECREATE");
	data->data_histograms[0]->Write();

	for (auto i = 0U; i < config.nparams; i++) {
		auto current_hist = pdfs->pdf_histograms[i];
		current_hist->Scale(results.popt[i]);
		current_hist->Write();
	}

	f->Close();
}

// @brief For now, simply plot the results (simple fit example)
auto ProcessResults(std::vector<NuFitData*> data, NuFitPDFs *pdfs, const NuFitConfig config, NuFitResults results) -> void {

}

}  // namespace NuFitter
