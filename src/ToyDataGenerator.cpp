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
// template <class Config>
auto generateToyData(const NuFitConfig config, const NuFitPDFs *pdfs) -> std::vector<NuFitData*> {

	std::vector<NuFitData*> data;

	for(auto t = 0U; t < config.ToyData; t++){

	std::vector<TH1D*> hists;
	std::vector<std::vector<double>> vec_data, bin_edges;
	std::vector<unsigned int> hist_ids;

	for (auto i = 1U; i <= config.data_hist_names.size(); i++) {
		if (std::find(config.hist_id.begin(), config.hist_id.end(), i) != config.hist_id.end()) {
				TH1D *hdata = new TH1D(config.data_hist_names[i-1].c_str(), config.data_hist_names[i-1].c_str(), config.nbins[i-1],
				    pdfs->bin_edges[i-1].front() + 1, pdfs->bin_edges[i-1].back());
				hists.push_back(hdata);
			}
		}

		gRandom = new TRandom3(0);
		gRandom->SetSeed(0);

		for (auto i = 0U; i < config.nSp_histos.size(); i++){
			for (auto k = 0U; k < config.nSp_histos[i]; k++) {
				auto idx = k + i * config.nSp_histos[0];
				auto current_hist = (TH1D*)pdfs->pdf_histograms[idx]->Clone();
				hists.at(i)->FillRandom(current_hist, config.param_initial_guess[idx]*config.param_eff[idx]);
			}
		}
		// Convert histogram to vector
			for (auto i = 1U; i <= config.data_hist_names.size(); i++) {
				std::vector<double> vec_data_hist;
				for (auto j = 0U; j <= config.nbins[i-1]; j++) {
					vec_data_hist.push_back(hists.at(i-1)->GetBinContent(j));
				}
				// Create bin_edges vector
				auto bin_edges_hist = getBinEdges_Toy(hists.at(i-1), config.nbins[i-1]);

				bin_edges.push_back(bin_edges_hist);
				vec_data.push_back(vec_data_hist);
				hist_ids.push_back(i-1);

			}


	auto *data_tofill = new NuFitData(vec_data, bin_edges, hists, hist_ids);
	data.push_back(data_tofill);

}
	return data;
}
}  // namespace NuFitter
