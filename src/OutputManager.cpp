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
#include "TTree.h"
// Project includes
#include "OutputManager.h"

struct Values {
	double mean, std_dev;
};

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

	// Create a std::vector<TH1D*> to fill with the fit results
    // TODO: bin width isn't always integer nor 1 !!
	auto range = config.emax - config.emin;

	std::vector<TH1D*> PDFsSum;
	for (auto i : data->hist_ids){
		auto name = "PDFsSum_" + config.data_hist_names[i];
		TH1D *PDFs_hists = new TH1D(name.c_str(), name.c_str(), range,
		    config.emin, config.emax);
		PDFsSum.push_back(PDFs_hists);
	}

	// Create a canvas to plot the results in
	auto nHists = data->hist_ids.size();
	TCanvas *c = new TCanvas("Results", "Results", 1500, 700);
	gPad->SetLogy();

	// Draw the data histograms
	TPad *padUp[nHists];
	TLegend *leg[nHists];
	for (auto i = 0U; i < nHists; i++){
		auto namePadUp = "PadUp_" + std::to_string(i+1);
		auto nameLeg = "Leg_" + std::to_string(i+1);
		padUp[i] = new TPad(namePadUp.c_str(), namePadUp.c_str(),
							i/static_cast<float>(nHists), 0.3,
							(1.+i)/static_cast<float>(nHists), 1.0);
		leg[i] = new TLegend(0.34,0.55 + i/5. ,0.54,0.85,NULL,"brNDC");
		c->cd();
		padUp[i]->Draw();
		padUp[i]->cd();
		gPad->SetLogy();
		data->data_histograms[i]->SetLineColor(kBlack);
		data->data_histograms[i]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
		data->data_histograms[i]->GetYaxis()->SetTitle("Events");
		data->data_histograms[i]->GetYaxis()->SetRangeUser(1, 1e5);
		data->data_histograms[i]->GetXaxis()->SetRangeUser(config.emin, config.emax);
		data->data_histograms[i]->Draw();
		leg[i]->SetTextAlign(13);
		leg[i]->SetTextSize(0.04);
		leg[i]->SetBorderSize(0);
		leg[i]->SetFillStyle(0);
		c->cd();
	}

	// Define colors to be used in the plots
	int *Colors = new int [13]{632,632,632,409,616,400,600,870,921,801,801,881,419};
	std::vector<int> colors;
	std::vector<TString> used_names;
	auto idx_col {0};

	// Plot the pdfs for each parameter
	for (auto i = 0U; i < results.paramVector.size(); i++) {
		auto parData = results.paramVector[i];
		for (auto el : parData) {
			auto j = el.idx_pdf;
			auto current_name = config.param_names[j];

			// Manage color used
			if (std::find(used_names.begin(), used_names.end(), current_name)
			    == used_names.end()) {
				used_names.push_back(current_name);
				colors.push_back(Colors[idx_col]);
			}

			auto current_hist = (TH1D*)pdfs->pdf_histograms[j]->Clone();
			current_hist->Scale(results.popt[i]/results.efficiencies[j]*config.param_eff[j]);

			// Update PDFSum
			for(auto k = 1U; k <= range; k++){
				PDFsSum.at(el.idx_hist-1)->SetBinContent(k,PDFsSum.at(el.idx_hist-1)->GetBinContent(k)+current_hist->GetBinContent(k-1+config.emin));
			}

			// Draw PDF
			padUp[el.idx_hist-1]->cd();
			current_hist->SetLineColor(colors[idx_col]);
			current_hist->SetMarkerColor(colors[idx_col]);
			current_hist->Draw("SAME");
			PDFsSum.at(el.idx_hist-1)->SetLineColor(632);
			PDFsSum.at(el.idx_hist-1)->SetMarkerColor(632);
			PDFsSum.at(el.idx_hist-1)->Draw("SAME");
			if(el.idx_hist == 2 && (current_name == "C11_2" || current_name == "C10" || current_name == "He6")){
				leg[el.idx_hist-1]->AddEntry(current_hist, config.param_names.at(j));
				leg[el.idx_hist-1]->Draw("SAME");
			}
			if(el.idx_hist == 1){
				leg[el.idx_hist-1]->AddEntry(current_hist, config.param_names.at(j));
				leg[el.idx_hist-1]->Draw("SAME");
			}
		}
		idx_col++;
		if (idx_col > 13) idx_col = 0;
	}

	// Fill vectors used to plot residuals
	std::vector<std::vector<double>> rec_energy;
	std::vector<std::vector<double>> residuals;
	for (auto i : data->hist_ids) {
		std::vector<double> res, rec;
		for(auto j = 0U; j < config.emax; j++){
			rec.push_back(j+config.emin);
			res.push_back((data->data_histograms[i]->GetBinContent(j+config.emin) -
			               PDFsSum[i]->GetBinContent(j+1)) /
						  sqrt(data->data_histograms[i]->GetBinContent(j+config.emin)));
		}
		rec_energy.push_back(rec);
		residuals.push_back(res);
	}

	// Plot residuals
	TPad *padDown[nHists];
	TGraph *Res[nHists];
	for (auto i = 0U; i < nHists; i++){
		auto namePadDown = "PadDown_" + std::to_string(i+1);
		padDown[i] = new TPad(namePadDown.c_str(), namePadDown.c_str(), i/static_cast<float>(nHists), 0, (1.+i)/static_cast<float>(nHists), 0.3);
		c->cd();
		padDown[i]->Draw();
		padDown[i]->cd();
		Res[i] = new TGraph(range, vec2Array(rec_energy[i]),
	                        vec2Array(residuals[i]));
		Res[i]->SetTitle("Residuals");
		Res[i]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
		Res[i]->GetYaxis()->SetTitle("(D-M)/sqrt(D)");
		Res[i]->GetYaxis()->CenterTitle(true);
		Res[i]->GetYaxis()->SetTitleSize(.05);
		Res[i]->GetXaxis()->SetTitleSize(.05);
		Res[i]->GetXaxis()->SetRangeUser(config.emin,config.emax);
		Res[i]->GetYaxis()->SetRangeUser(-4.,4.);
		Res[i]->SetLineWidth(1);
		Res[i]->Draw("AL");
	}
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
    // TODO: Same param on different hists can have different efficiencies!
    for (auto i = 0U; i < results.paramVector.size(); i++) {
        auto paramVec = results.paramVector[i];
        auto j = paramVec[0].idx_pdf;
        auto eff_exposure = factor / results.efficiencies[j];
        popt_cpd.push_back(results.popt[i] * eff_exposure);
        popt_err_cpd.push_back(popt_err[i] * eff_exposure);
    }

	// Write params with uncertainties in counts and cpd/kton
	outf << "------------\n" << "Fit results:\n" << "------------\n"
		 << "Species\tcounts\t\tsigma\t\trate(cpd/kton)\tsigma(cpd/kton)\n"
	     << std::scientific;
	outf.precision(4);
    for (auto i = 0U; i < results.paramVector.size(); i++) {  // For each parameter
        auto paramVec = results.paramVector[i];
        outf << config.param_names[paramVec[0].idx_pdf] << "\t";
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
	                const NuFitConfig config,
					std::vector<NuFitResults> results) -> void {

	auto factor {1. / (config.lifetime*config.mass_target)};


	//----------------------------------------------------------------
	//---------- Create the output rootfile --------------------------
	//----------------------------------------------------------------
	//----------------------------------------
	//----------- Plot the distributions -----
	//----------------------------------------
	auto root_filename = config.output_name + ".root";
	TFile *f = new TFile(root_filename.c_str(), "RECREATE");
	f->cd();

	TTree *tree = new TTree("Distributions","Distributions");

	std::vector<TString> names;
	int npar = results[0].paramVector.size();
	Values val[npar];

	//Loop to create all the branches of the TTree
	for (auto i = 0U; i < results[0].paramVector.size(); i++){
		auto paramVec = results[0].paramVector[i];
		auto name = config.param_names[paramVec[0].idx_pdf];
		names.push_back(name);
		tree->Branch(name, &val[i], "RecRates/D:StdDev/D");
	}

	//Loop to fill the TTree
	for (auto t = 0U; t < data.size(); t++){
		auto popt_err = results[t].getUncertainties();
		for (auto i = 0U; i < results[t].paramVector.size(); i++){
			auto paramVec = results[t].paramVector[i];
			auto j = paramVec[0].idx_pdf;
			auto eff_exposure = factor / results[t].efficiencies[j];
			val[i].mean = results[t].popt.at(i) * eff_exposure;
			val[i].std_dev = popt_err[i] * eff_exposure;
		}
		tree->Fill();
	}

	tree->Write();


	//------------------------------------------------------------
	//----------- Plot the results (only the last one) -----------
	//------------------------------------------------------------
	int lasttoy = data.size()-1;
	auto range = config.emax - config.emin;

	std::vector<TH1D*> PDFsSum;
	for (auto i : data[lasttoy]->hist_ids){
		auto name = "PDFsSum_" + config.data_hist_names[i];
		TH1D *PDFs_hists = new TH1D(name.c_str(), name.c_str(), range,
		config.emin, config.emax);
		PDFsSum.push_back(PDFs_hists);
	}

	// Create a canvas to plot the results in
	auto nHists = data[lasttoy]->hist_ids.size();
	TCanvas *c = new TCanvas("Plot", "Plot", 1500, 700);
	gPad->SetLogy();

	// Draw the data histograms
	TPad *padUp[nHists];
	TLegend *leg[nHists];
	for (auto i = 0U; i < nHists; i++){
		auto namePadUp = "PadUp_" + std::to_string(i+1);
		auto nameLeg = "Leg_" + std::to_string(i+1);
		padUp[i] = new TPad(namePadUp.c_str(), namePadUp.c_str(),
		                    i/static_cast<float>(nHists), 0.3,
							(1.+i)/static_cast<float>(nHists), 1.0);
		leg[i] = new TLegend(0.34,0.55 + i/5. ,0.54,0.85,NULL,"brNDC");
		c->cd();
		padUp[i]->Draw();
		padUp[i]->cd();
		gPad->SetLogy();
		data[lasttoy]->data_histograms[i]->SetLineColor(kBlack);
		data[lasttoy]->data_histograms[i]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
		data[lasttoy]->data_histograms[i]->GetYaxis()->SetTitle("Events");
		data[lasttoy]->data_histograms[i]->GetYaxis()->SetRangeUser(1, 1e5);
		data[lasttoy]->data_histograms[i]->GetXaxis()->SetRangeUser(config.emin, config.emax);
		data[lasttoy]->data_histograms[i]->Draw();
		leg[i]->SetTextAlign(13);
		leg[i]->SetTextSize(0.04);
		leg[i]->SetBorderSize(0);
		leg[i]->SetFillStyle(0);
		c->cd();
	}

	// Define colors to be used in the plots
	int *Colors = new int [13]{632,632,632,409,616,400,600,870,921,801,801,881,419};
	std::vector<int> colors;
	std::vector<TString> used_names;
	auto idx_col {0};

	// Plot the pdfs for each parameter
	for (auto i = 0U; i < results[lasttoy].paramVector.size(); i++) {
		auto parData = results[lasttoy].paramVector[i];
		for (auto el : parData) {
			auto j = el.idx_pdf;
			auto current_name = config.param_names[j];

			// Manage color used
			if (std::find(used_names.begin(), used_names.end(), current_name)
			    == used_names.end()) {
				used_names.push_back(current_name);
				colors.push_back(Colors[idx_col]);
			}

			auto current_hist = (TH1D*)pdfs->pdf_histograms[j]->Clone();
			current_hist->Scale(results[lasttoy].popt[i]/results[lasttoy].efficiencies[j]*config.param_eff[j]);

			// Update PDFSum
			for(auto k = 1U; k <= range; k++){
				PDFsSum.at(el.idx_hist-1)->SetBinContent(k,PDFsSum.at(el.idx_hist-1)->GetBinContent(k)+current_hist->GetBinContent(k-1+config.emin));
			}

			// Draw PDF
			padUp[el.idx_hist-1]->cd();
			current_hist->SetLineColor(colors[idx_col]);
			current_hist->SetMarkerColor(colors[idx_col]);
			current_hist->Draw("SAME");
			PDFsSum.at(el.idx_hist-1)->SetLineColor(632);
			PDFsSum.at(el.idx_hist-1)->SetMarkerColor(632);
			PDFsSum.at(el.idx_hist-1)->Draw("SAME");
			if(el.idx_hist == 2 && (current_name == "C11_2" || current_name == "C10" || current_name == "He6")){
				leg[el.idx_hist-1]->AddEntry(current_hist, config.param_names.at(j));
				leg[el.idx_hist-1]->Draw("SAME");
			}
			if(el.idx_hist == 1){
				leg[el.idx_hist-1]->AddEntry(current_hist, config.param_names.at(j));
				leg[el.idx_hist-1]->Draw("SAME");
			}
		}
		idx_col++;
		if (idx_col > 13) idx_col = 0;
	}


	// Fill vectors used to plot residuals
	std::vector<std::vector<double>> rec_energy;
	std::vector<std::vector<double>> residuals;
	for (auto i : data[lasttoy]->hist_ids) {
		std::vector<double> res, rec;
		for(auto j = 0U; j < config.emax; j++){
			rec.push_back(j+config.emin);
			res.push_back((data[lasttoy]->data_histograms[i]->GetBinContent(j+config.emin) -
			PDFsSum[i]->GetBinContent(j+1)) /
			sqrt(data[lasttoy]->data_histograms[i]->GetBinContent(j+config.emin)));
		}
		rec_energy.push_back(rec);
		residuals.push_back(res);
	}

	// Plot residuals
	TPad *padDown[nHists];
	TGraph *Res[nHists];

	for (auto i = 0U; i < nHists; i++){
		auto namePadDown = "PadDown_" + std::to_string(i+1);
		padDown[i] = new TPad(namePadDown.c_str(), namePadDown.c_str(), i/static_cast<float>(nHists), 0, (1.+i)/static_cast<float>(nHists), 0.3);
		c->cd();
		padDown[i]->Draw();
		padDown[i]->cd();
		Res[i] = new TGraph(range, vec2Array(rec_energy[i]),
		vec2Array(residuals[i]));
		Res[i]->SetTitle("Residuals");
		Res[i]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
		Res[i]->GetYaxis()->SetTitle("(D-M)/sqrt(D)");
		Res[i]->GetYaxis()->CenterTitle(true);
		Res[i]->GetYaxis()->SetTitleSize(.05);
		Res[i]->GetXaxis()->SetTitleSize(.05);
		Res[i]->GetXaxis()->SetRangeUser(config.emin,config.emax);
		Res[i]->GetYaxis()->SetRangeUser(-4.,4.);
		Res[i]->SetLineWidth(1);
		Res[i]->Draw("AL");
	}
	c->Write();
	f->Close();


	//----------------------------------------------------------------
	//---------- Create the output txt file --------------------------
	//----------------------------------------------------------------
	std::ofstream outf;
	auto out_filename = config.output_name + ".txt";
	outf.open(out_filename.c_str());

	for (auto t = 0U; t < data.size(); t++){

		outf << "------------------------------------------" << std::endl;
		outf << "Simulation number:\t" << t + 1 << std::endl;
		outf << "------------------------------------------" << std::endl;

		// Write the fit status
		outf << "[STATUS] Migrad status: ";
		if (results[t].errorflag != 0) {
			outf << "FAILED - ierflg =" << results[t].errorflag;
		} else {
			outf << "SUCCESS";
		}

		// Write the covariance matrix calculation status
		outf << std::endl << "[STATUS] Covariance matrix status: ";
		switch (results[t].errorflag_cov) {
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
			outf << "FAILED - istat=" << results[t].errorflag_cov;
		}
		outf << std::endl;

		// Convert counts to cpd/100t
		// cpd = count / (lifetime) / mass_target / efficiency;
		std::vector<double> popt_cpd, popt_err_cpd;
		auto popt_err = results[t].getUncertainties();
		// TODO: Same param on different hists can have different efficiencies!
		for (auto i = 0U; i < results[t].paramVector.size(); i++) {
			auto paramVec = results[t].paramVector[i];
			auto j = paramVec[0].idx_pdf;
			auto eff_exposure = factor / results[t].efficiencies[j];
			popt_cpd.push_back(results[t].popt[i] * eff_exposure);
			popt_err_cpd.push_back(popt_err[i] * eff_exposure);
		}

		// Write params with uncertainties in counts and cpd/kton
		outf << "------------\n" << "Fit results:\n" << "------------\n"
		<< "Species\tcounts\t\tsigma\t\trate(cpd/kton)\tsigma(cpd/kton)\n"
		<< std::scientific;
		outf.precision(4);
		for (auto i = 0U; i < results[t].paramVector.size(); i++) {  // For each parameter
			auto paramVec = results[t].paramVector[i];
			outf << config.param_names[paramVec[0].idx_pdf] << "\t";
			outf << results[t].popt[i] << "\t" << popt_err[i] << "\t";
			outf << popt_cpd[i] << "\t" << popt_err_cpd[i];
			outf << "\n";
		}

		// Write the covariance matrix
		outf << "------------------" << std::endl
		<< "Covariance matrix:" << std::endl
		<< "------------------" << std::endl;
		outf.precision(2);
		for (auto el : results[t].pcov) {
			for (auto sub : el) {
				outf << sub << "\t";
			}
			outf << std::endl;
		}
		// Write the correlation matrix
		auto test = results[t].getCorrMatrix();
		outf << "-------------------" << std::endl
		<< "Correlation matrix:" << std::endl
		<< "-------------------" << std::endl;
		for (auto el : test) {
			for (auto sub : el) {
				outf << sub << "\t";
			}
			outf << std::endl;
		}
	}
	// Write some information about the fit inputs
	outf << std::fixed
	<< "----------------------------------------------------------------------------------" << std::endl
	<< "Relevant fit inputs for the calculation:" << std::endl
	<< "----------------------------------------------------------------------------------" << std::endl
	<< "Exposure: " << config.lifetime << " days * "
	<< config.mass_target << " kton" << std::endl
	<< "Fit range: " << config.emin << "-" << config.emax << std::endl
	<< "Likelihood: '" << config.likelihood << "'" << std::endl;
	outf.close();
}

}  // namespace NuFitter
