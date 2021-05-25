//! @file	src/OutputManager.cpp
//! @brief	Implementation of class to process fit results
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard Includes
#include <iostream>
#include <fstream>
// ROOT includes
#include "TH1D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TCanvas.h"
// Project includes
#include "OutputManager.h"

namespace NuFitter {

// @brief For now, simply plot the results (simple fit example)
auto ProcessResults(NuFitData *data, NuFitPDFs *pdfs, const NuFitConfig config,
                    NuFitResults results) -> void {
	//----------------------------------------
	//----------- Plot the results -----------
	//----------------------------------------
	// Create a histogram with fit results
	// (maybe we can move it FitResults and make it become a member of it)
	TH1D* PDFsSum_histo = new TH1D("PDFsSum_histo","PDFsSum_histo", config.nbins, pdfs->bin_edges.front(), pdfs->bin_edges.back());

	// Open file to save the plots in
	auto root_filename = config.output_name + ".root";
	TFile *f = new TFile(root_filename.c_str(), "RECREATE");
	f->cd();

	TCanvas *c = new TCanvas("Results","Results",1500,700);
	gPad->SetLogy();

	TPad* Pad_up = new TPad("Pad_up", "Pad_up", 0, 0.3, 1.0, 1.0);
	Pad_up->Draw();
	Pad_up->cd();

	data->data_histograms[0]->SetLineColor(kBlack);
	data->data_histograms[0]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
	data->data_histograms[0]->GetYaxis()->SetTitle("Events");
	data->data_histograms[0]->Draw();

	//the maximum number of the species is 12 (Be7,pep,CNO,Bi210,K40,Kr85,U238,Th232,Po210,C10,He6,C11)
	int *Colors = new int [12]{632,632,409,616,400,600,870,921,632,801,881,419};

	TLegend *leg = new TLegend(0.54,0.55,0.74,0.85,NULL,"brNDC");
        leg->SetTextAlign(13);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

	for (auto i = 0U; i < config.npdfs; i++) {
		auto current_hist = pdfs->pdf_histograms[i];
		current_hist->SetLineColor(Colors[i]);
		current_hist->SetMarkerColor(Colors[i]);
		current_hist->Scale(results.popt[i]/results.efficiencies[i]);
		current_hist->Draw("SAME");

		for(auto j = 1U; j <= config.nbins; j++) {
			PDFsSum_histo->SetBinContent(j, PDFsSum_histo->GetBinContent(j)
			    +pdfs->pdf_histograms[i]->GetBinContent(j));
		}

		leg->AddEntry(current_hist,config.param_names.at(i));
		leg->Draw("SAME");
	}
	PDFsSum_histo->SetLineColor(632);
	PDFsSum_histo->SetMarkerColor(632);
	PDFsSum_histo->Draw("SAME");
	gPad->SetLogy();

	// Residuals
	double res[config.nbins];
	double rec_energy[config.nbins];

	for(auto i = 0U; i < config.nbins; i++){
		rec_energy[i] = i + pdfs->bin_edges.front();
		res[i] = (data->data_histograms[0]->GetBinContent(i)
		    - PDFsSum_histo->GetBinContent(i))
			/ sqrt(data->data_histograms[0]->GetBinContent(i));
	}

	c->cd();
  	TPad* Pad_down = new TPad("Pad_down", "Pad_down", 0.0, 0.0, 1.0, 0.3);
  	Pad_down->Draw();
  	Pad_down->cd();

	TGraph *Residuals = new TGraph(config.nbins,rec_energy,res);
	Residuals->SetTitle("Residuals");
	Residuals->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
	Residuals->GetYaxis()->SetTitle("(D-M)/sqrt(D)");
	Residuals->GetYaxis()->CenterTitle(true);
	Residuals->GetYaxis()->SetTitleSize(.05);
	Residuals->GetXaxis()->SetTitleSize(.05);
	Residuals->GetXaxis()->SetRangeUser(pdfs->bin_edges.front(),pdfs->bin_edges.back());
	Residuals->GetYaxis()->SetRangeUser(-4.,4.);
	Residuals->SetLineWidth(1);
	Residuals->Draw("AL");

	c->Write();
	f->Close();

	//----------------------------------------
	//------ Create the output txt file ------
	//----------------------------------------
	std::ofstream outf;
	auto out_filename = config.output_name + ".txt";
	outf.open(out_filename.c_str());

	// Write the fit status
	outf << "[STATUS] Migrad status: ";
	if (results.errorflag != 0) {
		outf << "FAILED - ierflg =" << results.errorflag;
	} else {
		outf << "SUCCESS";
	}

	// Write the covariance matrix calculation status
	outf << std::endl << "[STATUS] Covariance matrix status: ";
	switch (results.errorflag_cov) {
		case 0:
			outf << "FAILED - not calculated";
			break;
		case 1:
			outf << "FAILED - not accurate";
			break;
		case 2:
			outf << "FAILED - forced pos-def";
			break;
		case 3:
			outf << "SUCCESS";
			break;
		default:
			outf << "FAILED - istat=" << results.errorflag_cov;
	}
	outf << std::endl;

	// Convert counts to cpd/100t
	// cpd = count / (lifetime) / mass_target / efficiency;
	auto factor {1. / (config.lifetime*config.mass_target)};

	std::vector<double> popt_cpd, popt_err_cpd;
	auto popt_err = results.getUncertainties();
	for (auto i = 0U; i < config.npdfs; i++) {
		auto eff_exposure = factor / results.efficiencies[i];
		popt_cpd.push_back(results.popt[i] * eff_exposure);
		popt_err_cpd.push_back(popt_err[i] * eff_exposure);
	}

	// Write params with uncertainties in counts and cpd/kton
	outf << "Species\tcounts\t\tsigma\t\trate(cpd/kton)\tsigma(cpd/kton)\n"
	     << std::scientific;
	outf.precision(4);
	for (auto i = 0U; i < config.npdfs; i++) {
		outf << config.param_names[i] << "\t";
		outf << results.popt[i] << "\t" << popt_err[i] << "\t";
		outf << popt_cpd[i] << "\t" << popt_err_cpd[i];
		outf << "\n";
	}

	// Write the covariance matrix
	outf << "Covariance matrix:" << std::endl;
	outf.precision(2);
	for (auto el : results.pcov) {
		for (auto sub : el) {
			outf << sub << "\t";
		}
		outf << std::endl;
	}
	// Write the correlation matrix
	auto test = results.getCorrMatrix();
	outf << "Correlation matrix:" << std::endl;
	for (auto el : test) {
		for (auto sub : el) {
			outf << sub << "\t";
		}
		outf << std::endl;
	}
	outf.close();
}

// @brief Same as the other ProcessResults() but for toy data fit(s)
auto ProcessResults(std::vector<NuFitData*> data, NuFitPDFs *pdfs, const NuFitConfig config, NuFitResults results) -> void {

}

}  // namespace NuFitter
