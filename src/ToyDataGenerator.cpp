//! @file	src/ToyDataGenerator.cpp
//! @brief	Implementation of class to generate toy data
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard includes
#include <memory>
#include <iostream>
// ROOT includes
#include "TH1D.h"
#include "TRandom3.h"
// Project includes
#include "ToyDataGenerator.h"
#include "DataReader.h"
#include "TFile.h"

namespace NuFitter{

//Temporary. It's the same function as in DataReader. Maybe there is an easier way to useit without
//re-write it here
auto getBinEdges_Toy(TH1D *hist, unsigned int nbins) -> std::vector<double> {
	// Create bin_edges vector
	std::vector<double> bin_edges;
	for (auto i = 0U; i <= nbins; i++){
		bin_edges.push_back(hist->GetBinLowEdge(i));
	}
	bin_edges.push_back(hist->GetBinLowEdge(nbins) +
	                    hist->GetBinWidth(nbins));

	return bin_edges;
}

auto generateToyData(const NuFitConfig config, const NuFitPDFs *pdfs) -> std::vector<NuFitData*> {
	std::vector<NuFitData*> data;
	std::vector<TH1D*> histogr;
	auto nHists = config.nSp_histos_toy.size();

	// Fill the histogr vector with the TH1Ds in which the toy data will go
	for (auto i = 1U; i <= nHists; i++) {
		if (std::find(config.hist_id_toy.begin(), config.hist_id_toy.end(), i) != config.hist_id_toy.end()) {
			TH1D *hdata = new TH1D(config.data_hist_names[i-1].c_str(), config.data_hist_names[i-1].c_str(), config.nbins[i-1],
			pdfs->bin_edges[i-1].front() + 1, pdfs->bin_edges[i-1].back());
			histogr.push_back(hdata);
		}
	}

	// For each toy dataset
	for(auto t = 0U; t < config.ToyData; t++){
		// Make sure the TH1D(s) is(are) empty
		for (auto el : histogr) {
			el->Reset();
		}
		std::vector<std::vector<double>> vec_data, bin_edges;
		std::vector<unsigned int> hist_ids;
		gRandom = new TRandom3(0);
		gRandom->SetSeed(0);

		// Fill histogram from PDFs
		for (auto i = 0U; i < nHists; i++) {  // For each toy hist
			for (auto k = 0U; k < config.nSp_histos_toy[i]; k++) {  // For each pdf in hist
				auto idx = k + i * config.nSp_histos_toy[0];
				auto current_hist = (TH1D*)pdfs->pdf_histograms[idx]->Clone();
				auto n_expected = config.param_initial_guess_toy[idx]*config.param_eff_toy[idx];
				auto pois = gRandom->Poisson(n_expected);
				std::cout << n_expected << ", " << pois << std::endl;
				histogr.at(i)->FillRandom(current_hist, pois);
			}
		}

		// Convert histogram to vector
		for (auto i = 1U; i <= nHists; i++) {
			std::vector<double> vec_data_hist;
			for (auto j = 0U; j <= config.nbins[i-1]; j++) {
				vec_data_hist.push_back(histogr.at(i-1)->GetBinContent(j));
			}
			// Create bin_edges vector
			auto bin_edges_hist = getBinEdges_Toy(histogr.at(i-1), config.nbins[i-1]);

			bin_edges.push_back(bin_edges_hist);
			vec_data.push_back(vec_data_hist);
			hist_ids.push_back(i-1);
		}

		auto *data_tofill = new NuFitData(vec_data, bin_edges, histogr, hist_ids);
		data.push_back(data_tofill);

	}
	return data;
}
}  // namespace NuFitter
