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

template <class T>
auto vec2Array(std::vector<T> v) -> T* {
    T *array = &v[0];
    return array;
}

// @brief For now, simply plot the results (simple fit example)
auto ProcessResults(NuFitData *data, NuFitPDFs *pdfs, const NuFitConfig config,
                    NuFitResults results) -> void {
	//----------------------------------------
	//----------- Plot the results -----------
	//----------------------------------------
	// Open file to save the plots in
	auto root_filename = config.output_name + ".root";
	TFile *f = new TFile(root_filename.c_str(), "RECREATE");
	f->cd();

	TCanvas *c = new TCanvas("Results","Results",1500,700);
	gPad->SetLogy();

	//Histo_Sub
	c->cd();
	TPad* Pad_UpLeft = new TPad("Pad_UpLeft","Pad_UpLeft", 0., 0.3, 0.5, 1.0);
	Pad_UpLeft->Draw();
	Pad_UpLeft->cd();
	gPad->SetLogy();
	data->data_histograms[0]->SetLineColor(kBlack);
	data->data_histograms[0]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
	data->data_histograms[0]->GetYaxis()->SetTitle("Events");
	data->data_histograms[0]->Draw();

	//Histo_Tag
	// c->cd();
	// TPad* Pad_UpRight = new TPad("Pad_UpRight","Pad_UpRight", 0.5, 0.3, 1.0, 1.0);
	// Pad_UpRight->Draw();
	// Pad_UpRight->cd();
	// gPad->SetLogy();
	// data->data_histograms[1]->SetLineColor(kBlack);
	// data->data_histograms[1]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
	// data->data_histograms[1]->GetYaxis()->SetTitle("Events");
	// data->data_histograms[1]->Draw();

	//Legends
	TLegend *leg_UpLeft = new TLegend(0.34,0.55,0.54,0.85,NULL,"brNDC");
	leg_UpLeft->SetTextAlign(13);
	leg_UpLeft->SetTextSize(0.04);
	leg_UpLeft->SetBorderSize(0);
	leg_UpLeft->SetFillStyle(0);

	// TLegend *leg_UpRight = new TLegend(0.34,0.75,0.54,0.85,NULL,"brNDC");
	// leg_UpRight->SetTextAlign(13);
	// leg_UpRight->SetTextSize(0.04);
	// leg_UpRight->SetBorderSize(0);
	// leg_UpRight->SetFillStyle(0);

	//Be7,pep,Bi210,K40,Kr85,U238,Th232,Po210,C10,He6,C11)
	int *Colors = new int [11]{632,632,409,616,400,600,870,921,801,881,419};

	for (auto i = 0U; i < results.paramVector.size(); i++) {
		if (results.paramVector[i][0].idx_hist != 1) continue;
		auto j = results.paramVector[i][0].idx_pdf;

		auto current_hist = pdfs->pdf_histograms[j];
		current_hist->SetLineColor(Colors[j]);
		current_hist->SetMarkerColor(Colors[j]);
		current_hist->Scale(1./current_hist->Integral());
		current_hist->Scale(results.popt[i]/results.efficiencies[j]*config.param_eff[j]);
		Pad_UpLeft->cd();
		current_hist->Draw("SAME");
		leg_UpLeft->AddEntry(current_hist, config.param_names.at(j));
		leg_UpLeft->Draw("SAME");
	}


	// for (auto i : data->hist_ids) {
	// 	for(auto j = 0U; j < config.nSp_histos[i]; j++){
	//
	// 		if(i==0){
	// 			auto current_hist = pdfs->pdf_histograms[j];
	// 			current_hist->SetLineColor(Colors[j]);
	// 			current_hist->SetMarkerColor(Colors[j]);
	// 			current_hist->Scale(results.popt[j]/results.efficiencies[j]);
	// 			Pad_UpLeft->cd();
	// 			current_hist->Draw("SAME");
	// 			//PDFsSum[0]->Draw("SAME");
	// 			//gPad->SetLogy();
	// 			leg_UpLeft->AddEntry(current_hist,config.param_names.at(j));
	// 			leg_UpLeft->Draw("SAME");
	// 		}
	// 		// if(i==1){
	// 		// 	auto current_hist = pdfs->pdf_histograms[j+config.nSp_histos[0]];
	// 		// 	current_hist->SetLineColor(Colors[j]);
	// 		// 	current_hist->SetMarkerColor(Colors[j]);
	// 		// 	current_hist->Scale(results.popt[j]/results.efficiencies[j]);
	// 		// 	Pad_UpRight->cd();
	// 		// 	current_hist->Draw("SAME");
	// 		// 	//PDFsSum[1]->Draw("SAME");
	// 		// 	//gPad->SetLogy();
	// 		// 	Pad_UpRight->cd();
	// 		// 	leg_UpRight->AddEntry(current_hist,config.param_names.at(j+config.nSp_histos[0]));
	// 		// 	leg_UpRight->Draw("SAME");
	// 		// }
	// 	}
	// }


	// Residuals
	/*std::vector<std::vector<double>> residuals;
	std::vector<std::vector<double>> rec_energy;

	for (auto i : data->hist_ids) {
		for(auto j = 0U; j < config.nbins[i]; j++){
			rec_energy[i].push_back(j+pdfs->bin_edges[i].front());
			residuals[i].push_back((data->data_histograms[i]->GetBinContent(j)-PDFsSum[i]->GetBinContent(j))/sqrt(data->data_histograms[i]->GetBinContent(j)));
		}
	}

    c->cd();
    TPad* Pad_DownLeft = new TPad("Pad_DownLeft", "Pad_DownLeft", 0.0, 0.0, 0.5, 0.3);
    Pad_DownLeft->Draw();
    Pad_DownLeft->cd();

    // auto resArray = vec2Array(residuals[0]);
    TGraph *ResLeft = new TGraph(config.nbins[0], vec2Array(rec_energy[0]),
                                 vec2Array(residuals[0]));
    ResLeft->SetTitle("Residuals");
    ResLeft->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
    ResLeft->GetYaxis()->SetTitle("(D-M)/sqrt(D)");
    ResLeft->GetYaxis()->CenterTitle(true);
    ResLeft->GetYaxis()->SetTitleSize(.05);
    ResLeft->GetXaxis()->SetTitleSize(.05);
    ResLeft->GetXaxis()->SetRangeUser(pdfs->bin_edges[0].front(),pdfs->bin_edges[0].back());
    ResLeft->GetYaxis()->SetRangeUser(-4.,4.);
    ResLeft->SetLineWidth(1);
    ResLeft->Draw("AL");

    TPad* Pad_DownRight = new TPad("Pad_DownRight", "Pad_DownRight", 0.5, 0.0, 1.0, 0.3);
    Pad_DownRight->Draw();
    Pad_DownRight->cd();

    TGraph *ResRight = new TGraph(config.nbins[1], vec2Array(rec_energy[1]),
                                  vec2Array(residuals[1]));
    ResRight->SetTitle("Residuals");
    ResRight->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
    ResRight->GetYaxis()->SetTitle("(D-M)/sqrt(D)");
    ResRight->GetYaxis()->CenterTitle(true);
    ResRight->GetYaxis()->SetTitleSize(.05);
    ResRight->GetXaxis()->SetTitleSize(.05);
    ResRight->GetXaxis()->SetRangeUser(pdfs->bin_edges[1].front(),pdfs->bin_edges[1].back());
    ResRight->GetYaxis()->SetRangeUser(-4.,4.);
    ResRight->SetLineWidth(1);
    ResRight->Draw("AL");*/

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

	// Counts = PDF * parameter * param_eff
	// std::cout << "popt:" << std::endl;
	// for (auto i = 0U; i < results.paramVector.size(); i++) {
	// 	auto conv = results.popt[i]*config.param_eff[results.paramVector[i][0].idx_pdf]*factor/results.efficiencies[i];
	// 	std::cout << results.popt[i] << ", " << conv << ", " << config.param_eff[results.paramVector[i][0].idx_pdf] << ", " << results.efficiencies[i] << "; " << factor << std::endl;
	// }

	// Write params with uncertainties in counts and cpd/kton
	outf << "------------\n" << "Fit results:\n" << "------------\n"
		 << "Species\tcounts\t\tsigma\t\trate(cpd/kton)\tsigma(cpd/kton)\n"
	     << std::scientific;
	outf.precision(4);
	for (auto i = 0U; i < config.npdfs; i++) {
		outf << config.param_names[i] << "\t";
		outf << results.popt[i] << "\t" << popt_err[i] << "\t";
		outf << popt_cpd[i] << "\t" << popt_err_cpd[i];
		outf << "\n";
	}

	// Write the covariance matrix
	outf << "------------------" << std::endl
	     << "Covariance matrix:" << std::endl
	     << "------------------" << std::endl;
	outf.precision(2);
	for (auto el : results.pcov) {
		for (auto sub : el) {
			outf << sub << "\t";
		}
		outf << std::endl;
	}
	// Write the correlation matrix
	auto test = results.getCorrMatrix();
	outf << "-------------------" << std::endl
	     << "Correlation matrix:" << std::endl
	     << "-------------------" << std::endl;
	for (auto el : test) {
		for (auto sub : el) {
			outf << sub << "\t";
		}
		outf << std::endl;
	}

	// Write some information about the fit inputs
	outf << std::fixed
	     << "----------------------------------------" << std::endl
	     << "Relevant fit inputs for the calculation:" << std::endl
		 << "----------------------------------------" << std::endl
	     << "Exposure: " << config.lifetime << " days * "
		 << config.mass_target << " kton" << std::endl
		 << "Fit range: " << config.emin << "-" << config.emax << std::endl
		 << "Likelihood: '" << config.likelihood << "'" << std::endl;


	outf.close();
}

// @brief Same as the other ProcessResults() but for toy data fit(s)
auto ProcessResults(std::vector<NuFitData*> data, NuFitPDFs *pdfs,
	                const NuFitConfig config, NuFitResults results) -> void {

}

}  // namespace NuFitter
