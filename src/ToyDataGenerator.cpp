//! @file	src/ToyDataGenerator.cpp
//! @brief	Implementation of class to generate toy data
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard includes
#include <memory>
#include <iostream>
// ROOT includes
#include "TH1D.h"
// Project includes
#include "ToyDataGenerator.h"
#include "DataReader.h"
#include "TFile.h"

#define INT_MAX 2147483647

namespace NuFitter {

//Temporary. It's the same function as in DataReader. Maybe there is an easier way to useit without
//re-write it here
auto getBinEdges_(TH1D *hist, unsigned int nbins) -> std::vector<double> {
	// Create bin_edges vector
	std::vector<double> bin_edges;
	for (auto i = 0U; i <= nbins; i++){
		bin_edges.push_back(hist->GetBinLowEdge(i));
	}
	bin_edges.push_back(hist->GetBinLowEdge(nbins) +
	                    hist->GetBinWidth(nbins));

	return bin_edges;
}

// @brief Constructor for NuFitToyData
// pass reference to pdf pointer -> no need for delete later
NuFitToyData::NuFitToyData(NuFitPDFs *&pdfs_) {
	pdfs = pdfs_;

	// Fill the histogr vector with the TH1Ds in which the toy data will go
	auto nHists = config->nSp_histos_toy.size();
	for (auto i = 1U; i <= nHists; i++) {
		if (std::find(config->hist_id_toy.begin(), config->hist_id_toy.end(), i)
		    != config->hist_id_toy.end()) {

			hdata = new TH1D(config->data_hist_names[i-1].c_str(),
			                 config->data_hist_names[i-1].c_str(), config->nbins[i-1],
			pdfs->bin_edges[i-1].front() + 1, pdfs->bin_edges[i-1].back());
			histogr.push_back(hdata);
			// Get and save the bin edges for re-use later
			bin_edges.push_back(getBinEdges_(histogr.at(i-1), config->nbins[i-1]));
			hist_ids.push_back(i-1);
		}
	}
}

// @brief Destructor for NuFitToyData
// @brief No need to delete pdfs -> pointer was passed by reference
NuFitToyData::~NuFitToyData() {
	if (dataset) delete dataset;
	// This also implicitely does "delete hdata" that points to histogr.back()
	for (auto el : histogr) {
		delete el;
	}
}

// @brief Create dataset on the fly by filling the histograms with scaled PDFs
auto NuFitToyData::loadDataset(unsigned int idx_dataset) -> void {
	// Initialise
	std::vector<std::vector<double>> vec_data;
	assert(idx_dataset < config->ToyData);

	// Make sure the TH1D(s) is(are) empty
	for (auto el : histogr) {
		el->Reset();
	}

	// Fill histogram from PDFs
	auto samples = config->param_sampled[idx_dataset];
	for (auto parData : config->paramVector_toy) {  // For each parameter
		for (auto el : parData) {  // For each PDF
			auto j = el.idx_pdf;
			auto current_hist = (TH1D*)pdfs->pdf_histograms[j]->Clone();
			auto n_samples = samples[j];
			// Take care of problems when n_samples is greater than
			// the numeric limit for integers
			if (n_samples > INT_MAX) std::cout << "Warning: the rate of " << config->param_names[j] << " is very high, it could cause numeric problems." << std::endl;
			while (n_samples > INT_MAX) {
				histogr[el.idx_hist-1]->FillRandom(current_hist, INT_MAX);
				n_samples -= INT_MAX;
			}
			histogr[el.idx_hist-1]->FillRandom(current_hist, n_samples);
			delete current_hist;
		}
	}

	// TODO: add a check that no bin contains more than DOUBLE_MAX?

	// Convert histogram to vector
	// Loop of data histograms only if they were found in the toy specieslist!
	for (auto i = 1U; i <= config->data_hist_names.size(); i++) {
		if (std::find(config->hist_id_toy.begin(), config->hist_id_toy.end(), i)
		    != config->hist_id_toy.end()) {

			std::vector<double> vec_data_hist;
			for (auto j = 0U; j <= config->nbins[i-1]; j++) {
				vec_data_hist.push_back(histogr.at(i-1)->GetBinContent(j));
			}
			vec_data.push_back(vec_data_hist);
		}
	}

	if (dataset) delete dataset;
	// Only the last one needs root histograms to be plotted
	if (idx_dataset < config->ToyData-1) {
		dataset = new NuFitData(vec_data, bin_edges, hist_ids);
	} else {
		dataset = new NuFitData(vec_data, bin_edges, histogr, hist_ids);
	}

}
}  // namespace NuFitter
