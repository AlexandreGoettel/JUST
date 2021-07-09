//! @file	src/DataReader.cpp
//! @brief	Implementation of classes to read data and PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard includes
#include <memory>
#include <iostream>
// Project includes
#include "DataReader.h"
#include "Parser.h"


namespace NuFitter {

// TODO: do we want this function somewhere else?
auto getBinEdges(TH1D *hist, unsigned int nbins) -> std::vector<double> {
	// Create bin_edges vector
	std::vector<double> bin_edges;
	for (auto i = 0U; i <= nbins; i++){
		bin_edges.push_back(hist->GetBinLowEdge(i));
	}
	bin_edges.push_back(hist->GetBinLowEdge(nbins) +
	                    hist->GetBinWidth(nbins));

	return bin_edges;
}

// @brief Constructor for NuFitPDFs
NuFitPDFs::NuFitPDFs() {
	std::cout << "Constructing NuFitPDFs" << std::endl;
}

// @brief Destructor for NuFitPDFs
NuFitPDFs::~NuFitPDFs() {
	std::cout << "Destructing NuFitPDFs" << std::endl;
	if (file_pdf) file_pdf->Close();
	// Histograms should already be deleted by root when file_pdf is closed
	// for (auto el : pdf_histograms) {
	// 	delete el;
	// }
}

// @brief function to read pdfs from files
auto NuFitPDFs::Read(const std::vector<TString>& pdf_names, const unsigned int& npdfs,
	                 const std::vector<unsigned int>& hist_id) -> void {
	// Read the PDFs
	std::cout << "INFO: Reading PDFs from " << config->pdf_name << std::endl;
	file_pdf = new TFile(config->pdf_name.c_str());
	for (auto i = 0U; i < npdfs; i++) {
		auto pdf = (TH1D*)file_pdf->Get(pdf_names[i]);
		// Make sure PDFs are normalised to 1
		pdf->Scale(1./pdf->Integral());  // TODO: include bin width?
		pdf_histograms.push_back(pdf);
	}

	// Convert histograms to vectors
	for (auto n = 0U; n < npdfs; n++) {
		std::vector<double> current_pdf;
		auto ref_pdf = pdf_histograms[n];
		for (auto i = 0U; i <= config->nbins[hist_id[n]-1]; i++) {
			current_pdf.push_back(ref_pdf->GetBinContent(i));
		}
		pdfs.push_back(current_pdf);
	}

	// Create bin_edges vector for each used data histogram
	// Only if the histogram was found in specieslist
	for (auto i = 1U; i <= config->data_hist_names.size(); i++) {
		if (std::find(hist_id.begin(), hist_id.end(), i) != hist_id.end()) {
			// TODO: different histograms can have different binning!
			auto bin_edges_tmp = getBinEdges(pdf_histograms[0], config->nbins[i-1]);
			bin_edges.push_back(bin_edges_tmp);
		}
	}
}

// @brief Constructor for NuFitData including TH1Ds for plotting
NuFitData::NuFitData(std::vector<std::vector<double>> data_,
	                 std::vector<std::vector<double>> bin_edges_,
	                 std::vector<TH1D*> data_histograms_,
                     std::vector<unsigned int> ids_) {
	data = data_;
	bin_edges = bin_edges_;
	data_histograms = data_histograms_;
	hist_ids = ids_;
}

// @brief Constructor for NuFitData not including TH1Ds
NuFitData::NuFitData(std::vector<std::vector<double>> data_,
	                 std::vector<std::vector<double>> bin_edges_,
                     std::vector<unsigned int> ids_) {
	data = data_;
	bin_edges = bin_edges_;
	hist_ids = ids_;
}

namespace Data {

auto Read() -> NuFitData* {
	// Initialise variables
	std::vector<TH1D*> hists;
	std::vector<std::vector<double>> vec_data, bin_edges;
	std::vector<unsigned int> hist_ids;

	// Read the histograms if their id was found in specieslist
	std::cout << "INFO: Reading data from " << config->data_name << std::endl;
	for (auto i = 1U; i <= config->data_hist_names.size(); i++) {
		if (std::find(config->hist_id.begin(), config->hist_id.end(), i) != config->hist_id.end()) {
			// Read histogram
			TFile *file_data = new TFile(config->data_name.c_str());
			TH1D* hdata = (TH1D*)file_data->Get(config->data_hist_names[i-1].c_str());

			// Convert histogram to vector
			std::vector<double> vec_data_hist;
			for (auto j = 0U; j <= config->nbins[i-1]; j++) {
				vec_data_hist.push_back(hdata->GetBinContent(j));
			}

			// Create bin_edges vector
			auto bin_edges_hist = getBinEdges(hdata, config->nbins[i-1]);

			// Save to vectors
			hists.push_back(hdata);
			bin_edges.push_back(bin_edges_hist);
			vec_data.push_back(vec_data_hist);
			hist_ids.push_back(i-1);
	    }
	}

	std::cout << "[DATA READER]" << "nHists: " << vec_data.size()
	          << ", nBins: " << config->nbins[0] << std::endl;

	// Create and return NuFitData object
	auto *data = new NuFitData(vec_data, bin_edges, hists, hist_ids);
	return data;
}

}  // namespace Data
}  // namespace NuFitter
