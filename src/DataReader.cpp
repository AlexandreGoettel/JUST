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

	TFile *file_data = new TFile(config.data_name.c_str());
	TH1D* hdata = (TH1D*)file_data->Get(config.histo_data.c_str());

	// Convert histogram to vector
	std::vector<double> vec_data;
	for (auto i = 1U; i <= config.nbins; i++) {
		vec_data.push_back(hdata->GetBinContent(i));
	}

	// Create bin_edges vector
	auto bin_edges = getBinEdges(hdata, config.nbins);

	// Also save histograms (for plotting later)
	std::vector<TH1D*> hists;
	hists.push_back(hdata);

	// Create and return NuFitData object
	auto *data = new NuFitData(vec_data, bin_edges, hists);
	return data;
}

}  // namespace Data

namespace PDFs {

auto Read(const NuFitConfig config) -> NuFitPDFs* {
	// Generate the PDFs
		TFile *file_pdf = new TFile(config.pdf_name.c_str());
		TH1D** hPDF = new TH1D*[config.npdfs];

		hPDF[config.npdfs-config.npdfs] = (TH1D*)file_pdf->Get(config.histo_pdf_signal.c_str());
		//hPDF[config.npdfs-config.npdfs]->Scale(1./hPDF[config.npdfs-config.npdfs]->Integral());
		hPDF[config.npdfs-config.npdfs+1] = (TH1D*)file_pdf->Get(config.histo_pdf_background.c_str());
		//hPDF[config.npdfs-config.npdfs+1]->Scale(1./hPDF[config.npdfs-config.npdfs+1]->Integral());


	// Convert histograms to vector
	std::vector<double> vec_background, vec_signal;
	for (auto i = 1U; i <= config.nbins; i++) {
		vec_background.push_back(hPDF[config.npdfs-config.npdfs+1]->GetBinContent(i));
		vec_signal.push_back(hPDF[config.npdfs-config.npdfs]->GetBinContent(i));
	}
	std::vector<std::vector<double>> pdfs;
	pdfs.push_back(vec_signal);
	pdfs.push_back(vec_background);

	// Create bin_edges vector
	auto bin_edges = getBinEdges(hPDF[config.npdfs-config.npdfs+1], config.nbins);

	// Also save histograms in vector (for plotting later)
	std::vector<TH1D*> hists;
	hists.push_back(hPDF[config.npdfs-config.npdfs]);
	hists.push_back(hPDF[config.npdfs-config.npdfs+1]);

	// Create and return NuFitPDFs object with the variables
	auto *output = new NuFitPDFs(pdfs, bin_edges, hists);
	return output;
}

}  // namespace PDFs
}  // namespace NuFitter
