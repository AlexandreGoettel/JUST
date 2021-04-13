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
// Project includes
#include "DataReader.h"
#include "Parser.h"


namespace NuFitter {

template <class T, class C, class ... params>
auto genHist(T function, C functionName, C histName,
	         C histDescription, int nBins, double xmin, double xmax,
			 int nEvents, params ... p) -> TH1D* {
	// Create a TF1 object representing the function
	TF1 *f = new TF1(functionName, function, xmin, xmax);
	f->SetParameters(p ...);

	// Randomize seed
	gRandom = new TRandom();
	gRandom->SetSeed(0);

	// Create and fill a histogram with that function
	TH1D *h = new TH1D(histName, histDescription, nBins, xmin, xmax);
	h->FillRandom(functionName, nEvents);

	return h;
}

NuFitPDFs::NuFitPDFs(std::vector<std::vector<double>> pdfs_,
	                 std::vector<double> bin_edges_,
				     std::vector<TH1D*> pdf_histograms_) {
	pdfs = pdfs_;
	bin_edges = bin_edges_;
	pdf_histograms = pdf_histograms_;
}

NuFitData::NuFitData(std::vector<double>data_, std::vector<double> bin_edges_,
	                 std::vector<TH1D*> data_histograms_) {
	data = data_;
	bin_edges = bin_edges_;
	data_histograms = data_histograms_;
}

auto getBinEdges(TH1D *hist, unsigned int nbins) -> std::vector<double> {
	// Create bin_edges vector
	std::vector<double> bin_edges;
	for (auto i = 1U; i <= nbins; i++){
		bin_edges.push_back(hist->GetBinLowEdge(i));
	}
	bin_edges.push_back(hist->GetBinLowEdge(nbins) +
	                    hist->GetBinWidth(nbins));

	return bin_edges;
}

namespace Data {

auto Read(const NuFitConfig config) -> NuFitData* {
	// TODO
	// ##################################
	// This is all temporary for testing!
	// ##################################
	// Generate the data
	auto nEvents_data {15000};
	auto h_total = genHist("expo(0) + .5*gaus(2)", "f_tot", "h_tot",
	                       "Data", config.nbins, 0, 10, nEvents_data,
						   0.2, -1, 1, 2, .5);

	// Convert histogram to vector
	std::vector<double> vec_data;
	for (auto i = 1U; i <= config.nbins; i++) {
		vec_data.push_back(h_total->GetBinContent(i));
	}

	// Create bin_edges vector
	auto bin_edges = getBinEdges(h_total, config.nbins);

	// Also save histograms (for plotting later)
	std::vector<TH1D*> hists;
	hists.push_back(h_total);

	// Create and return NuFitData object
	auto *data = new NuFitData(vec_data, bin_edges, hists);
	return data;
}

}  // namespace Data

namespace PDFs {

auto Read(const NuFitConfig config) -> NuFitPDFs* {
	// TODO
	// ##################################
	// This is all temporary for testing!
	// ##################################
	// Generate the PDFs
	auto nEvents_pdfs {100000};
	auto h_background = genHist("expo(0)", "f_background", "h_background",
                       "Exponential background", config.nbins, 0, 10,
                       nEvents_pdfs, 0.2, -1);
   	h_background->Scale(1./nEvents_pdfs);
   	auto h_signal = genHist("gaus(0)", "f_signal", "h_signal",
   	                   "Gaussian signal", config.nbins, 0, 10, nEvents_pdfs,
   				       1, 2, .5);
   	h_signal->Scale(1./nEvents_pdfs);

	// Convert histograms to vector
	std::vector<double> vec_background, vec_signal;
	for (auto i = 1U; i <= config.nbins; i++) {
		vec_background.push_back(h_background->GetBinContent(i));
		vec_signal.push_back(h_signal->GetBinContent(i));
	}
	std::vector<std::vector<double>> pdfs;
	pdfs.push_back(vec_signal);
	pdfs.push_back(vec_background);

	// Create bin_edges vector
	auto bin_edges = getBinEdges(h_background, config.nbins);

	// Also save histograms in vector (for plotting later)
	std::vector<TH1D*> hists;
	hists.push_back(h_signal);
	hists.push_back(h_background);

	// Create and return NuFitPDFs object with the variables
	auto *output = new NuFitPDFs(pdfs, bin_edges, hists);
	return output;
}

}  // namespace PDFs
}  // namespace NuFitter
