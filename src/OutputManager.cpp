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
#include "Parser.h"


namespace NuFitter {

// @brief output values to be inserted in a ROOT tree branch
struct Values {
	double fit_rate, fit_rate_err, fit_counts_tot, fit_counts_range;
};
struct ValuesToy {
	double toy_rate, gen_counts;
};
struct ValuesParam {
	double injected_rate, initial_guess, stepsize, lowerlim, upperlim;
	unsigned int isFixed;
};

// @brief Convert a vector to an array, works because C++ guarantees ordering
template <class T>
auto vec2Array(std::vector<T> v) -> T* {
    T *array = &v[0];
    return array;
}

// // @brief OutputManager constructor for data fits
// OutputManager::OutputManager(NuFitData*& data_, NuFitPDFs*& pdfs_, NuFitResults*& results_) {
// 	data = data_;
// 	pdfs = pdfs_;
// 	std::vector<NuFitResults*> res {results_};
// 	results = res;
// 	// TODO??
// }
//
// // @brief OutputManager constructor for toy fits
// OutputManager::OutputManager(NuFitToyData*& data_, NuFitPDFs*& pdfs_,
// 	                         NuFitPDFs*& pdfs_toy_, std::vector<NuFitResults*>& results_) {
// 	data_toy = data_;
// 	pdfs = pdfs_;
// 	pdfs_toy = pdfs_toy_;
// 	results = results_;
// }

// @brief Convert a results parameter vector of counts into cpd/kton
auto toCpdPerkton(std::vector<double> ref, NuFitResults *&results) -> std::vector<double> {
	assert(ref.size() == results->paramVector.size());
	auto factor {1. / config->exposure};
	std::vector<double> output;

	for (auto i = 0U; i < results->paramVector.size(); i++) {
		// TODO: Same param on different hists can have different efficiencies!
		auto j = results->paramVector[i][0].idx_pdf;
		auto eff_exposure = factor;
		if (config->param_fixed[j] != 1) {
			// Fixed parameters are already scaled to the full spectrum
			eff_exposure /= results->efficiencies[j];
		}
		output.push_back(ref[i] * eff_exposure);
	}
	return output;
}

// @brief Destructor for OutputManager
OutputManager::~OutputManager() {
	if (f->IsOpen()) f->Close();
	delete PDFsSum;
}

// @brief Open the root file to be used to write the results in
auto OutputManager::initRootFile(std::string &filename) -> void {
	f = new TFile(filename.c_str(), "RECREATE");
	f->cd();
}

// @brief Close the root file
auto OutputManager::closeRootFile() -> void {
	f->Close();
}

// @brief Initialise and fill PDFsSum so it can be shared later
auto OutputManager::makePDFsSum(NuFitData *&data, NuFitPDFs*& pdfs,
                                NuFitResults *& results) -> void {
	// Temporary
	// TODO: bin width isn't always integer nor 1 !!
	auto range = config->emax - config->emin;

	PDFsSum = new std::vector<TH1D*>;
	for (auto i : data->hist_ids){
		auto name = "PDFsSum_" + config->data_hist_names[i];
		TH1D *PDFs_hists = new TH1D(name.c_str(), name.c_str(), range,
		    config->emin, config->emax);
		PDFsSum->push_back(PDFs_hists);
	}

	// Fill PDFsSum
	for (auto i = 0U; i < results->paramVector.size(); i++) {
		auto parData = results->paramVector[i];
		for (auto el : parData) {
			auto j = el.idx_pdf;
			auto n = el.idx_hist-1;
			auto current_name = config->param_names[j];
			auto current_hist = (TH1D*)pdfs->pdf_histograms[j]->Clone();
			auto current_PDFsSum = PDFsSum->at(n);

			if (config->param_fixed[j] == 1) {
				current_hist->Scale(results->popt[i] * config->param_eff[j]);
			} else {
				current_hist->Scale(results->popt[i] * config->param_eff[j]
				                    / results->efficiencies[j]);
			}

			// Update PDFSum
			for(auto k = 1U; k <= current_PDFsSum->GetNbinsX(); k++){
				current_PDFsSum->SetBinContent(k, current_PDFsSum->GetBinContent(k)
				    + current_hist->GetBinContent(k - 1 + config->emin));
			}
		}
	}

	paramVector = results->paramVector;
}

// @brief Create, fill, and save a tree with some meta information
auto OutputManager::writeParamTree(NuFitResults *&results) -> void {
	f->cd();
	f->mkdir("parameters");
	f->cd("parameters");
	TTree *paramTree = new TTree("Parameters", "Parameters");
	TTree *configTree = new TTree("Config", "Config");

	// Initialise branches
	auto npar = results->paramVector.size();
	ValuesParam val[npar];
	for (auto i = 0U; i < npar; i++) {
		// Initialise
		auto paramVec = results->paramVector[i];
		auto j = paramVec[0].idx_pdf;
		auto name = config->param_names[j];
		paramTree->Branch(name, &val[i], "injected_rate/D:initial_guess/D:stepsize/D:lowerlim/D:upperlim/D:isFixed/i");

		// Check if its a toy or data fit
		if (config->param_initial_guess_toy.size() > 0) {
			val[i].injected_rate = config->param_initial_guess_toy[j]/config->exposure;
		}

		// Fill
		// double injected_rate, initial_guess, stepsize, lowerlim, upperlim
		val[i].initial_guess = config->param_initial_guess[j]/config->exposure;
		val[i].stepsize = config->param_stepsize[j]/config->exposure;
		val[i].lowerlim = config->param_lowerlim[j]/config->exposure;
		val[i].upperlim = config->param_upperlim[j]/config->exposure;
		val[i].isFixed = config->param_fixed[j];
	}
	paramTree->Fill();

	// Now init tree for non-parameter numbers
	double exposure, emin, emax;
	unsigned int seed;
	configTree->Branch("exposure", &exposure, "exposure/D");
	configTree->Branch("emin", &emin, "emin/D");
	configTree->Branch("emax", &emax, "emax/D");
	configTree->Branch("seed", &seed, "emax/i");

	// Fill
	exposure = config->exposure;
	emin = config->emin;
	emax = config->emax;
	seed = config->seed;
	configTree->Fill();

	// Write to file
	paramTree->Write();
	configTree->Write();
	f->cd();
}

// @brief Create, fill, and write a tree containing the fit results
auto OutputManager::writeFitTree(std::vector<NuFitResults*>& results) -> void {
	// Initialise
	TTree *fitTree = new TTree("FitResults", "FitResults");
	auto npar = results[0]->paramVector.size();
	Values val[npar];

	// Create all the branches of the Tree
	for (auto i = 0U; i < results[0]->paramVector.size(); i++) {
		auto paramVec = results[0]->paramVector[i];
		auto name = config->param_names[paramVec[0].idx_pdf];
		fitTree->Branch(name, &val[i], "fit_rate/D:fit_rate_err/D:fit_counts_tot/D:fit_counts_range/D");
	}

	// Fill the Tree
	for (auto t = 0U; t < results.size(); t++) {  // For each toy dataset
		auto results_ = results[t];
		auto popt_cpd = toCpdPerkton(results_->popt, results_);
		auto popt_err_cpd = toCpdPerkton(results_->getUncertainties(), results_);

		// Get the fit result information
		for (auto i = 0U; i < npar; i++) {
			val[i].fit_rate = popt_cpd[i];
			val[i].fit_rate_err = popt_err_cpd[i];
			val[i].fit_counts_range = results_->popt[i];
			val[i].fit_counts_tot = popt_cpd[i]*config->exposure;
		}

		fitTree->Fill();
	}
	f->cd();
	fitTree->Write();
}

// @brief Create, fill, and write a tree containing toy generation information
auto OutputManager::writeToyTree() -> void {
	// Initialise
	TTree *toyTree = new TTree("ToyGeneration", "ToyGeneration");

	// Create branches for the toy data tree
	auto nparToy = config->paramVector_toy.size();
	ValuesToy valToy[nparToy];
	for (auto i = 0U; i < nparToy; i++) {
		auto j = config->paramVector_toy[i][0].idx_pdf;
		auto name = config->pdf_names_toy[j];
		toyTree->Branch(name, &valToy[i], "toy_rate/D:gen_counts/D");
	}

	// Fill the TTrees
	for (auto t = 0U; t < config->ToyData; t++) {  // For each toy dataset
		// Get the toy data information
		auto samples = config->param_sampled[t];
		for (auto i = 0U; i < config->paramVector_toy.size(); i++) {  // For each param
			auto gen_counts {0UL};
			auto parData = config->paramVector_toy[i];

			// Number of PDFs controlled by the parameter
			auto nPDFs_param = parData.size();
			for (auto j = 0U; j < nPDFs_param; j++) {
				auto k = parData[j].idx_pdf;
				// Dividing by param_eff points back directly to the total
				// number of physical events. So take the average over the
				// pdfs to get the toy-generated poisson-fluctuated rate
				gen_counts += samples[k] / config->param_eff[k] / nPDFs_param;
			}

			valToy[i].toy_rate = gen_counts/config->exposure;
			valToy[i].gen_counts = gen_counts;
		}
		toyTree->Fill();
	}
	f->cd();
	toyTree->Write();
}

// @brief Write the output of one fit to file
auto OutputManager::fitToFile(std::ofstream &outf, NuFitResults *&results) -> void {
	// Write the fit status
	outf << "[STATUS] Migrad status: ";
	if (results->errorflag != 0) {
		outf << "FAILED - ierflg =" << results->errorflag;
	} else {
		outf << "SUCCESS";
	}

	// Write the covariance matrix calculation status
	outf << std::endl << "[STATUS] Covariance matrix status: ";
	switch (results->errorflag_cov) {
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
			outf << "FAILED - istat=" << results->errorflag_cov;
	}
	outf << std::endl;

	// Convert counts to cpd/100t
	// cpd = count / (lifetime) / mass_target / efficiency;
	auto popt_err = results->getUncertainties();
	auto popt_cpd = toCpdPerkton(results->popt, results);
	auto popt_err_cpd = toCpdPerkton(popt_err, results);

	// Write params with uncertainties in counts and cpd/kton
	outf << "------------\n" << "Fit results:\n" << "------------\n"
		 << "Species\tcounts\t\tsigma\t\trate(cpd/kton)\tsigma(cpd/kton)\n"
	     << std::scientific;
	outf.precision(4);
    for (auto i = 0U; i < results->paramVector.size(); i++) {  // For each parameter
        auto paramVec = results->paramVector[i];
        outf << config->param_names[paramVec[0].idx_pdf] << "\t";
        outf << results->popt[i] << "\t" << popt_err[i] << "\t";
        outf << popt_cpd[i] << "\t" << popt_err_cpd[i];
		outf << "\n";
    }

	// Write the covariance matrix
	outf << "------------------" << std::endl
	     << "Covariance matrix:" << std::endl
	     << "------------------" << std::endl;
	outf.precision(2);
	for (auto el : results->pcov) {
		for (auto sub : el) {
			outf << sub << "\t";
		}
		outf << std::endl;
	}
	// Write the correlation matrix
	auto test = results->getCorrMatrix();
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
	     << "Exposure: " << config->lifetime << " days * "
		 << config->mass_target << " kton" << std::endl
		 << "Fit range: " << config->emin << "-" << config->emax << std::endl
		 << "Likelihood: '" << config->likelihood << "'" << std::endl
		 << "Random seed: " << config->seed << std::endl;
}

auto OutputManager::plotToFile(NuFitData *&data, NuFitPDFs *&pdfs,
	                           NuFitResults *&results) -> void {
	f->cd();
	// Alternative: build PDFsSum with OututManager constructors?
	if (!PDFsSum) makePDFsSum(data, pdfs, results);

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
		data->data_histograms[i]->GetXaxis()->SetRangeUser(config->emin, config->emax);
		data->data_histograms[i]->Draw();
		leg[i]->SetTextAlign(13);
		leg[i]->SetTextSize(0.04);
		leg[i]->SetBorderSize(0);
		leg[i]->SetFillStyle(0);
		c->cd();
	}

	// Define colors to be used in the pdfs
	int *Colors = new int [13]{632,632,632,409,616,400,600,870,921,801,801,881,419};
	auto idx_col {0};

	// Plot the pdfs for each parameter
	for (auto i = 0U; i < results->paramVector.size(); i++) {
		auto parData = results->paramVector[i];
		for (auto el : parData) {
			auto j = el.idx_pdf;
			auto n = el.idx_hist-1;
			auto current_name = config->param_names[j];
			auto current_hist = (TH1D*)pdfs->pdf_histograms[j]->Clone();
			auto current_PDFsSum = PDFsSum->at(n);

			if (config->param_fixed[j] == 1) {
				current_hist->Scale(results->popt[i] * config->param_eff[j]);
			} else {
				current_hist->Scale(results->popt[i] * config->param_eff[j]
				                    / results->efficiencies[j]);
			}

			// Draw PDF
			padUp[n]->cd();
			current_hist->SetLineColor(Colors[idx_col]);
			current_hist->SetMarkerColor(Colors[idx_col]);
			current_hist->Draw("SAME");
			current_PDFsSum->SetLineColor(632);
			current_PDFsSum->SetMarkerColor(632);
			current_PDFsSum->Draw("SAME");
			if(n == 1 && (current_name == "C11_2" || current_name == "C10" || current_name == "He6")){
				leg[n]->AddEntry(current_hist, config->param_names.at(j));
				leg[n]->Draw("SAME");
			}
			if(n == 0){
				leg[n]->AddEntry(current_hist, config->param_names.at(j));
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
		for(auto j = 0U; j < config->emax; j++){
			rec.push_back(j+config->emin);
			res.push_back((data->data_histograms[i]->GetBinContent(j+config->emin) -
			               PDFsSum->at(i)->GetBinContent(j+1)) /
						  sqrt(data->data_histograms[i]->GetBinContent(j+config->emin)));
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
		Res[i] = new TGraph(PDFsSum->at(i)->GetNbinsX(), vec2Array(rec_energy[i]),
	                        vec2Array(residuals[i]));
		Res[i]->SetTitle("Residuals");
		Res[i]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
		Res[i]->GetYaxis()->SetTitle("(D-M)/sqrt(D)");
		Res[i]->GetYaxis()->CenterTitle(true);
		Res[i]->GetYaxis()->SetTitleSize(.05);
		Res[i]->GetXaxis()->SetTitleSize(.05);
		Res[i]->GetXaxis()->SetRangeUser(config->emin,config->emax);
		Res[i]->GetYaxis()->SetRangeUser(-4.,4.);
		Res[i]->SetLineWidth(1);
		Res[i]->Draw("AL");
	}
	c->Write();
}

// @brief Function to draw extra PDFs for comparison
auto OutputManager::drawPDFs(NuFitPDFs *&pdfs_fit, NuFitPDFs *&pdfs_toy) -> void {
	f->cd();

	TCanvas *c = new TCanvas("PDFs", "PDFs", 1500, 700);
	gPad->SetLogy();
	TPad *pad[2];
	// TODO: auto-adjust to the number of histograms..
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

	for (auto i = 0U; i < paramVector.size(); i++) {
		auto parData = paramVector[i];
		for (auto el : parData) {
			auto j = el.idx_pdf;
			auto current_name = config->param_names[j];
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

	// Plot the pdfs_toy for each parameter
	auto idx_col_toy {0};
	for (auto i = 0U; i < paramVector.size(); i++) {
		auto parData = paramVector[i];
		for (auto el : parData) {
			auto j = el.idx_pdf;
			auto current_name = config->param_names[j];
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
auto ProcessResults(NuFitData *&data, NuFitPDFs *&pdfs,
	                NuFitResults *&results) -> void {
	auto manager = OutputManager();

	//----------------------------------------
	//------ Create the output rootfile ------
	//----------------------------------------
	auto root_filename = config->output_name + ".root";
	manager.initRootFile(root_filename);

	//----------------------------------------
	//------ Create and fill the TTrees ------
	//----------------------------------------
	manager.writeParamTree(results);
	std::vector<NuFitResults*> results_vec {results};
	manager.writeFitTree(results_vec);

	//----------------------------
	//---------- Plots -----------
	//----------------------------
	// Save example plot
	manager.plotToFile(data, pdfs, results);
	manager.closeRootFile();

	//----------------------------------------
	//------ Create the output txt file ------
	//----------------------------------------
	std::ofstream outf;
	auto out_filename = config->output_name + ".txt";
	outf.open(out_filename.c_str());

	manager.fitToFile(outf, results);
	outf.close();
}

// @brief Same as the other ProcessResults() but for toy data fit(s)
auto ProcessResults(NuFitToyData *&data, NuFitPDFs *&pdfs_toy, NuFitPDFs *&pdfs,
					std::vector<NuFitResults*> &results) -> void {
	auto manager = OutputManager();

	//----------------------------------------
	//------ Create the output rootfile ------
	//----------------------------------------
	auto root_filename = config->output_name + ".root";
	manager.initRootFile(root_filename);

	//----------------------------------------
	//------ Create and fill the TTrees ------
	//----------------------------------------
	manager.writeParamTree(results[0]);
	manager.writeFitTree(results);
	manager.writeToyTree();

	//----------------------------
	//---------- Plots -----------
	//----------------------------
	// Save example plot
	manager.plotToFile(data->dataset, pdfs, results.back());
	// Draw pdfs separatly for comparison purposes (once per file only)
	manager.drawPDFs(pdfs, pdfs_toy);

	manager.closeRootFile();

	//----------------------------------------
	//------ Create the output txt file ------
	//----------------------------------------
	std::ofstream outf;
	auto out_filename = config->output_name + ".txt";
	outf.open(out_filename.c_str());

	for (auto t = 0U; t < config->ToyData; t++){  // For each toy dataset
		outf << "------------------------------------------" << std::endl;
		outf << "Simulation number:\t" << t + 1 << std::endl;
		outf << "------------------------------------------" << std::endl;
		manager.fitToFile(outf, results[t]);
	}
	outf.close();
}

}  // namespace NuFitter
