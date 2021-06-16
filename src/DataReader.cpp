//! @file	src/DataReader.cpp
//! @brief	Implementation of classes to read data and PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard includes
#include <memory>
#include <iostream>
// ROOT includes
#include "TF1.h"
#include "TRandom3.h"
#include <TFile.h>
// Project includes
#include "DataReader.h"
#include "Parser.h"


namespace NuFitter {
NuFitPDFs::NuFitPDFs(std::vector<std::vector<double>> pdfs_,
	                 std::vector<std::vector<double>> bin_edges_,
				     std::vector<TH1D*> pdf_histograms_) {
	pdfs = pdfs_;
	bin_edges = bin_edges_;
	pdf_histograms = pdf_histograms_;
}

NuFitData::NuFitData(std::vector<std::vector<double>> data_,
	                 std::vector<std::vector<double>> bin_edges_,
	                 std::vector<TH1D*> data_histograms_,
                     std::vector<unsigned int> ids_) {
	data = data_;
	bin_edges = bin_edges_;
	data_histograms = data_histograms_;
	hist_ids = ids_;
}

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

namespace Data {

auto Read(const NuFitConfig config) -> NuFitData* {
	// Initialise variables
	std::vector<TH1D*> hists;
	std::vector<std::vector<double>> vec_data, bin_edges;
	std::vector<unsigned int> hist_ids;

	// Read the histograms if their id was found in specieslist
	std::cout << "INFO: Reading data from " << config.data_name << std::endl;
	for (auto i = 1U; i <= config.data_hist_names.size(); i++) {
		if (std::find(config.hist_id.begin(), config.hist_id.end(), i) != config.hist_id.end()) {
			// Read histogram
			TFile *file_data = new TFile(config.data_name.c_str());
			TH1D* hdata = (TH1D*)file_data->Get(config.data_hist_names[i-1].c_str());

			// Convert histogram to vector
			std::vector<double> vec_data_hist;
			for (auto j = 0U; j <= config.nbins[i-1]; j++) {
				vec_data_hist.push_back(hdata->GetBinContent(j));
			}

			// Create bin_edges vector
			auto bin_edges_hist = getBinEdges(hdata, config.nbins[i-1]);

			// Save to vectors
			hists.push_back(hdata);
			bin_edges.push_back(bin_edges_hist);
			vec_data.push_back(vec_data_hist);
			hist_ids.push_back(i-1);
	    }
	}
	
	std::cout << "[DATA READER]" << "nHists: " << vec_data.size()
	          << ", nBins: " << config.nbins[0] << std::endl;

	// Create and return NuFitData object
	auto *data = new NuFitData(vec_data, bin_edges, hists, hist_ids);
	return data;
}

}  // namespace Data

namespace PDFs {

auto Read(const NuFitConfig config) -> NuFitPDFs* {
	// Read the PDFs and convert them to vectors
	std::cout << "INFO: Reading PDFs from " << config.pdf_name << std::endl;
	TFile *file_pdf = new TFile(config.pdf_name.c_str());
	std::vector<TH1D*> hPDFs;
	for (auto i = 0U; i < config.npdfs; i++) {
		hPDFs.push_back((TH1D*)file_pdf->Get(config.pdf_names[i]));
	}

	// Convert histograms to vectors
	std::vector<std::vector<double>> pdfs;
	for (auto n = 0U; n < config.npdfs; n++) {
		std::vector<double> current_pdf;
		for (auto i = 0U; i <= config.nbins[config.hist_id[n]-1]; i++) {
			current_pdf.push_back(hPDFs[n]->GetBinContent(i));
		}
		pdfs.push_back(current_pdf);
	}

	// Create bin_edges vector for each used data histogram
	// Only if the histogram was found in specieslist
	std::vector<std::vector<double>> bin_edges;
	for (auto i = 1U; i <= config.data_hist_names.size(); i++) {
		if (std::find(config.hist_id.begin(), config.hist_id.end(), i) != config.hist_id.end()) {
			// TODO: different histograms can have different binning!
			auto bin_edges_tmp = getBinEdges(hPDFs[0], config.nbins[i-1]);
			bin_edges.push_back(bin_edges_tmp);
		}
	}

	std::cout << "[PDF READER]" << "nPDF: " << pdfs.size() << std::endl
	          << "[PDF READER]" << "nHists: " << bin_edges.size() << std::endl;

	// Create and return NuFitPDFs object with the variables
	auto *output = new NuFitPDFs(pdfs, bin_edges, hPDFs);
	return output;
}

}  // namespace PDFs
}  // namespace NuFitter
