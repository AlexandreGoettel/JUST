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

	//Create a std::vector<TH1D*> to fill with the fit results
	std::vector<TH1D*> PDFsSum;

	for (auto i : data->hist_ids){
		auto name = "PDFsSum_" + config.data_hist_names[i];
		TH1D *PDFs_hists = new TH1D(name.c_str(), name.c_str(), config.nbins[i], pdfs->bin_edges[i].front(), pdfs->bin_edges[i].back());
		PDFsSum.push_back(PDFs_hists);
	}

	TCanvas *c = new TCanvas("Results","Results",1500,700);
	gPad->SetLogy();

	TPad *padUp[data->hist_ids.size()];
	TLegend *leg[data->hist_ids.size()];

	for (int i = 0; i < data->hist_ids.size(); i++){
		auto namePadUp = "PadUp_" + std::to_string(i+1);
		auto nameLeg = "Leg_" + std::to_string(i+1);
		padUp[i] = new TPad(namePadUp.c_str(), namePadUp.c_str(), 0. + i/2., 0.3, 0.5 + i/2., 1.0);
		leg[i] = new TLegend(0.34,0.55 + i/5. ,0.54,0.85,NULL,"brNDC");
		c->cd();
		padUp[i]->Draw();
		padUp[i]->cd();
		gPad->SetLogy();
		data->data_histograms[i]->SetLineColor(kBlack);
		data->data_histograms[i]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
		data->data_histograms[i]->GetYaxis()->SetTitle("Events");
		data->data_histograms[i]->GetYaxis()->SetRangeUser(1,1e5);
		data->data_histograms[i]->Draw();
		leg[i]->SetTextAlign(13);
		leg[i]->SetTextSize(0.04);
		leg[i]->SetBorderSize(0);
		leg[i]->SetFillStyle(0);
		c->cd();
	}

		//Be7,pep,Bi210,Kr85,Po210,U238,Th232,K40,C11,C11_2,C10,He6
		int *Colors = new int [12]{632,632,409,616,400,600,870,921,801,801,881,419};
		std::vector<int> colors;
		std::vector<TString> used_names;
		auto idx_col {0};

    for (auto i = 0U; i < results.paramVector.size(); i++) {
        auto parData = results.paramVector[i];
        for (auto el : parData) {
            auto j = el.idx_pdf;
						auto current_name = config.param_names[el.idx_pdf];
						if (std::find(used_names.begin(), used_names.end(), current_name)
							== used_names.end()) {
								used_names.push_back(current_name);
								colors.push_back(Colors[idx_col]);
							}

            auto current_hist = (TH1D*)pdfs->pdf_histograms[j]->Clone();
            current_hist->Scale(results.popt[i]/results.efficiencies[j]*config.param_eff[j]);

          if (el.idx_hist == 1) {
						for(auto k = 1U; k <= config.nbins[el.idx_hist-1]; k++){
						PDFsSum.at(el.idx_hist-1)->SetBinContent(k,PDFsSum.at(el.idx_hist-1)->GetBinContent(k)+current_hist->GetBinContent(k));
						}

    				padUp[el.idx_hist-1]->cd();
						current_hist->SetLineColor(colors[idx_col]);
            current_hist->SetMarkerColor(colors[idx_col]);
    				current_hist->Draw("SAME");
						PDFsSum.at(el.idx_hist-1)->SetLineColor(632);
						PDFsSum.at(el.idx_hist-1)->SetMarkerColor(632);
						PDFsSum.at(el.idx_hist-1)->Draw("SAME");
    				leg[el.idx_hist-1]->AddEntry(current_hist, config.param_names.at(j));
    				leg[el.idx_hist-1]->Draw("SAME");

    			}
					if (el.idx_hist == 2) {
						for(auto k = 1U; k <= config.nbins[el.idx_hist-1]; k++){
							PDFsSum.at(el.idx_hist-1)->SetBinContent(k,PDFsSum.at(el.idx_hist-1)->GetBinContent(k)+current_hist->GetBinContent(k));
						}

          	padUp[el.idx_hist-1]->cd();
						current_hist->SetLineColor(colors[idx_col]);
            current_hist->SetMarkerColor(colors[idx_col]);
          	current_hist->Draw("SAME");
						PDFsSum.at(el.idx_hist-1)->SetLineColor(632);
						PDFsSum.at(el.idx_hist-1)->SetMarkerColor(632);
						PDFsSum.at(el.idx_hist-1)->Draw("SAME");
						if(j==17||j==18||j==19){

    				leg[el.idx_hist-1]->AddEntry(current_hist, config.param_names.at(j));
    				leg[el.idx_hist-1]->Draw("SAME");
					}
        }
      }
			idx_col++;
    }

	// Residuals
	std::vector<std::vector<double>> rec_energy;
	std::vector<std::vector<double>> residuals;

	for (auto i : data->hist_ids) {
		std::vector<double> res,rec;
		for(auto j = 0U; j < config.nbins[i]; j++){
				rec.push_back(j+pdfs->bin_edges[i].front());
				res.push_back((data->data_histograms[i]->GetBinContent(j)-PDFsSum[i]->GetBinContent(j))/sqrt(data->data_histograms[i]->GetBinContent(j)));
		}
		rec_energy.push_back(rec);
		residuals.push_back(res);
	}

	TPad *padDown[data->hist_ids.size()];
	TGraph *Res[data->hist_ids.size()];

	for (int i = 0; i < data->hist_ids.size(); i++){
		auto namePadDown = "PadDown_" + std::to_string(i+1);
		padDown[i] = new TPad(namePadDown.c_str(),namePadDown.c_str(), 0. + i/2., 0., 0.5 + i/2., 0.3);
		c->cd();
		padDown[i]->Draw();
		padDown[i]->cd();
		Res[i] = new TGraph(config.nbins[i], vec2Array(rec_energy[i]),
	                                 vec2Array(residuals[i]));
		Res[i]->SetTitle("Residuals");
		Res[i]->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
		Res[i]->GetYaxis()->SetTitle("(D-M)/sqrt(D)");
		Res[i]->GetYaxis()->CenterTitle(true);
		Res[i]->GetYaxis()->SetTitleSize(.05);
		Res[i]->GetXaxis()->SetTitleSize(.05);
		Res[i]->GetXaxis()->SetRangeUser(pdfs->bin_edges[i].front(),pdfs->bin_edges[i].back());
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
	                const NuFitConfig config, NuFitResults results) -> void {

}

}  // namespace NuFitter
