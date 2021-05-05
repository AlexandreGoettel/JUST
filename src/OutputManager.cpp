//! @file	src/OutputManager.cpp
//! @brief	Implementation of class to process fit results
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Includes
#include <iostream>
#include "TH1D.h"
#include "TFile.h"
#include "TStyle.h"
#include <TLegend.h>
#include <TGraph.h>
#include "TCanvas.h"
#include "OutputManager.h"

namespace NuFitter {

// @brief For now, simply plot the results (simple fit example)
auto ProcessResults(NuFitData *data, NuFitPDFs *pdfs, const NuFitConfig config, NuFitResults results) -> void {

	// Create a histogram with fit results (maybe we can move it FitResults and make it become a member of it)
	TH1D* PDFsSum_histo = new TH1D("PDFsSum_histo","PDFsSum_histo", config.nbins, pdfs->bin_edges.front(), pdfs->bin_edges.back());

	// Open file to save the plots in
	TFile *f = new TFile(config.output_name.c_str(), "RECREATE");
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

	for (auto i = 0U; i < config.nparams; i++) {
		auto current_hist = pdfs->pdf_histograms[i];
		current_hist->SetLineColor(Colors[i]);
		current_hist->SetMarkerColor(Colors[i]);
		current_hist->Scale(results.popt[i]);
		current_hist->Draw("SAME");

		for(auto j = 1U; j <= config.nbins; j++)	PDFsSum_histo->SetBinContent(j,PDFsSum_histo->GetBinContent(j)+pdfs->pdf_histograms[i]->GetBinContent(j));

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
		res[i] = (data->data_histograms[0]->GetBinContent(i)-PDFsSum_histo->GetBinContent(i))/sqrt(data->data_histograms[0]->GetBinContent(i));
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
}

// @brief For now, simply plot the results (simple fit example)
auto ProcessResults(std::vector<NuFitData*> data, NuFitPDFs *pdfs, const NuFitConfig config, NuFitResults results) -> void {

}

}  // namespace NuFitter
