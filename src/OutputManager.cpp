//! @file	src/OutputManager.cpp
//! @brief	Implementation of class to process fit results
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard Includes
#include <iostream>
// ROOT includes
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
// Project includes
#include "OutputManager.h"


namespace NuFitter {

// @brief output values to be inserted in a ROOT tree branch
struct Values {
	double inj_rate, fit_rate, fit_rate_err, fit_counts_tot, fit_counts_range;
};
struct ValuesToy {
	double inj_rate_toy, toy_rate;
	unsigned int gen_counts;
	//double toy_rate;
};

template <class T>
auto vec2Array(std::vector<T> v) -> T* {
    T *array = &v[0];
    return array;
}

// @brief Convert a results parameter vector of counts into cpd/kton
auto toCpdPerkton(std::vector<double> ref, const NuFitConfig config,
	NuFitResults results) -> std::vector<double> {
	assert(ref.size() == results.paramVector.size());
	auto factor {1. / config.exposure};
	std::vector<double> output;

	// TODO: Same param on different hists can have different efficiencies!
	for (auto i = 0U; i < results.paramVector.size(); i++) {
		auto paramVec = results.paramVector[i];
		auto j = paramVec[0].idx_pdf;
		auto eff_exposure = factor;
		if (config.param_fixed[j] != 1) {
			// Fixed parameters are already scaled to the full spectrum
			eff_exposure /= results.efficiencies[j];
		}
		output.push_back(ref[i] * eff_exposure);
	}
	return output;
}

// @brief Write the output of one fit to file
auto fitToFile(std::ofstream &outf, NuFitResults results,
	                const NuFitConfig config) -> void {
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
	auto popt_err = results.getUncertainties();
	auto popt_cpd = toCpdPerkton(results.popt, config, results);
	auto popt_err_cpd = toCpdPerkton(popt_err, config, results);

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
}

auto plotToFile(TFile *f,  NuFitData *data, NuFitPDFs *pdfs,
	            const NuFitConfig config, NuFitResults results) -> void {
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
		data->data_histograms[i]->SetLineColor(kBlack);
		data->data_histograms[i]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
		data->data_histograms[i]->GetYaxis()->SetTitle("Events");
		// TODO: can't use hard-coded values
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
	auto idx_col {0};

	// Plot the pdfs for each parameter
	for (auto i = 0U; i < results.paramVector.size(); i++) {
		auto parData = results.paramVector[i];
		for (auto el : parData) {
			auto j = el.idx_pdf;
			auto n = el.idx_hist-1;
			auto current_name = config.param_names[j];

			auto current_hist = (TH1D*)pdfs->pdf_histograms[j]->Clone();
			if (config.param_fixed[j] == 1) {
				current_hist->Scale(results.popt[i] * config.param_eff[j]);
			} else {
				current_hist->Scale(results.popt[i] * config.param_eff[j]
				                    / results.efficiencies[j]);
			}

			// Update PDFSum
			for(auto k = 1U; k <= range; k++){
				PDFsSum.at(n)->SetBinContent(k,PDFsSum.at(n)->GetBinContent(k)+current_hist->GetBinContent(k-1+config.emin));
			}

			// Draw PDF
			padUp[n]->cd();
			current_hist->SetLineColor(Colors[idx_col]);
			current_hist->SetMarkerColor(Colors[idx_col]);
			current_hist->Draw("SAME");
			PDFsSum.at(n)->SetLineColor(632);
			PDFsSum.at(n)->SetMarkerColor(632);
			PDFsSum.at(n)->Draw("SAME");
			if(n == 1 && (current_name == "C11_2" || current_name == "C10" || current_name == "He6")){
				leg[n]->AddEntry(current_hist, config.param_names.at(j));
				leg[n]->Draw("SAME");
			}
			if(n == 0){
				leg[n]->AddEntry(current_hist, config.param_names.at(j));
				leg[n]->Draw("SAME");
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
}

auto DrawPDFs(TFile *f, NuFitPDFs *pdfs_toy, NuFitPDFs *pdfs_fit, const NuFitConfig config, NuFitResults results) -> void {

		f->cd();

		TCanvas *c = new TCanvas("PDFs", "PDFs", 1500, 700);
		gPad->SetLogy();
		TPad *pad[2];
		pad[0] = new TPad("PDFs_model", "PDFs_model",0.,0.,0.5,1.0);
		pad[1] = new TPad("PDFs_data", "PDFs_data",0.5,0.,1.0,1.0);
		pad[0]->Draw();
		gPad->SetLogy();
		pad[1]->Draw();
		gPad->SetLogy();
		c->cd();

		// Define colors to be used in the plots
		int *Colors = new int [13]{632,632,632,409,616,400,600,870,921,801,801,881,419};
		auto idx_col {0};

		// Plot the pdfs_fit for each parameter
		for (auto i = 0U; i < results.paramVector.size(); i++) {
			auto parData = results.paramVector[i];
			for (auto el : parData) {
				auto j = el.idx_pdf;
				auto n = el.idx_hist-1;
				auto current_name = config.param_names[j];

				auto current_hist = (TH1D*)pdfs_fit->pdf_histograms[j]->Clone();

				// Draw PDF
				pad[0]->cd();
				gPad->SetLogy();
				current_hist->SetLineColor(Colors[idx_col]);
				current_hist->SetMarkerColor(Colors[idx_col]);
				current_hist->SetTitle("PDFs: model");
				current_hist->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
				current_hist->GetXaxis()->SetRangeUser(0,10000);
				current_hist->GetYaxis()->SetRangeUser(1,0.02);
				current_hist->Draw("SAME");

			}
			idx_col++;
			if (idx_col > 13) idx_col = 0;
		}

		auto idx_col_toy {0};

		for (auto i = 0U; i < results.paramVector.size(); i++) {
			auto parData = results.paramVector[i];
			for (auto el : parData) {
				auto j = el.idx_pdf;
				auto n = el.idx_hist-1;
				auto current_name = config.param_names[j];

				auto current_hist = (TH1D*)pdfs_toy->pdf_histograms[j]->Clone();

				// Draw PDF
				pad[1]->cd();
				gPad->SetLogy();
				current_hist->SetLineColor(Colors[idx_col_toy]);
				current_hist->SetMarkerColor(Colors[idx_col_toy]);
				current_hist->SetTitle("PDFs: data");
				current_hist->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
				current_hist->GetYaxis()->SetTitle("Events");
				current_hist->GetXaxis()->SetRangeUser(0,10000);
				current_hist->GetYaxis()->SetRangeUser(1,0.02);
				current_hist->Draw("SAME");

			}
			idx_col_toy++;
			if (idx_col_toy > 13) idx_col_toy = 0;
		}
		c->Write();

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

	plotToFile(f, data, pdfs, config, results);
	f->Close();

	//----------------------------------------
	//------ Create the output txt file ------
	//----------------------------------------
	std::ofstream outf;
	auto out_filename = config.output_name + ".txt";
	outf.open(out_filename.c_str());

	fitToFile(outf, results, config);
	outf.close();
}

// @brief Same as the other ProcessResults() but for toy data fit(s)
auto ProcessResults(NuFitToyData* data, NuFitPDFs *pdfs_toy, NuFitPDFs *pdfs,
	                const NuFitConfig config,
					std::vector<NuFitResults> results) -> void {
	//----------------------------------------
	//------ Create the output rootfile ------
	//----------------------------------------
	auto root_filename = config.output_name + ".root";
	TFile *f = new TFile(root_filename.c_str(), "RECREATE");
	f->cd();

	DrawPDFs(f, pdfs_toy, pdfs, config, results[0]);

	// Initialise TTrees
	TTree *fitTree = new TTree("Distributions", "Distributions");
	TTree *toyTree = new TTree("ToyGeneration", "ToyGeneration");
	auto npar = results[0].paramVector.size();
	Values val[npar];

	// Create all the branches of the TTrees
	// fit_rate, fit_rate_err, fit_counts_tot, fit_counts_range, gen_counts
	for (auto i = 0U; i < results[0].paramVector.size(); i++){
		auto paramVec = results[0].paramVector[i];
		auto name = config.param_names[paramVec[0].idx_pdf];
		fitTree->Branch(name, &val[i], "inj_rate/D:fit_rate/D:fit_rate_err/D:fit_counts_tot/D:fit_counts_range/D");
	}

	// Create branches for the toy data tree
	auto nparToy = config.paramVector_toy.size();
	ValuesToy valToy[nparToy];
	for (auto i = 0U; i < nparToy; i++) {
		auto j = config.paramVector_toy[i][0].idx_pdf;
		auto name = config.pdf_names_toy[j];
		toyTree->Branch(name, &valToy[i], "inj_rate_toy/D:toy_rate/D:gen_counts/i");
	}

	// Fill the TTrees
	for (auto t = 0U; t < config.ToyData; t++){  // For each toy dataset
		auto results_ = results[t];
		auto popt_cpd = toCpdPerkton(results_.popt, config, results_);
		auto popt_err_cpd = toCpdPerkton(results_.getUncertainties(), config, results_);

		// Get the fit result information
		for (auto i = 0U; i < npar; i++){
			val[i].inj_rate = config.param_initial_guess[i]/config.exposure;
			val[i].fit_rate = popt_cpd[i];
			val[i].fit_rate_err = popt_err_cpd[i];
			val[i].fit_counts_range = results_.popt[i];
			val[i].fit_counts_tot = popt_cpd[i]*config.exposure;
		}

		// Get the toy data information
		auto samples = config.param_sampled[t];
		for (auto i = 0U; i < config.paramVector_toy.size(); i++) {
			auto gen_counts {0U};
			auto parData = config.paramVector_toy[i];
			for (auto el : parData) {
				gen_counts += samples[el.idx_pdf];
			}

			valToy[i].inj_rate_toy = config.param_initial_guess_toy[i]/config.exposure;
			valToy[i].toy_rate = gen_counts/config.exposure;
			valToy[i].gen_counts = gen_counts;
		}
		fitTree->Fill();
		toyTree->Fill();
	}
	fitTree->Write();
	toyTree->Write();

	//----------------------------------------
	//---------- Save example plot -----------
	//----------------------------------------
	plotToFile(f, data->dataset, pdfs, config, results.back());
	f->Close();

	//----------------------------------------
	//------ Create the output txt file ------
	//----------------------------------------
	std::ofstream outf;
	auto out_filename = config.output_name + ".txt";
	outf.open(out_filename.c_str());

	for (auto t = 0U; t < config.ToyData; t++){  // For each toy dataset
		outf << "------------------------------------------" << std::endl;
		outf << "Simulation number:\t" << t + 1 << std::endl;
		outf << "------------------------------------------" << std::endl;
		fitToFile(outf, results[t], config);
	}
	outf.close();
}

}  // namespace NuFitter
