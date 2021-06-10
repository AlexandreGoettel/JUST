//! @file	src/Parser.cpp
//! @brief	Implementation of classes to parse cmdl args and config files
//! @author	Alexandre Göttel
//! @date	2021-03-18

//============================================================================
// Includes
// Standard includes
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>

// ROOT includes
#include "TFile.h"
#include "TH1D.h"
// Project includes
#include "Parser.h"

//============================================================================
// Implementations

namespace NuFitter{

// @brief problems in opening or reading input files
auto ErrorReading(std::string& filename) -> void {
	std::ifstream setfile(filename);
	if(setfile.fail()){
		std::cout << "Opening " << filename << " for reading.\n";
		std::cout << "The "<< filename <<" file cannot be opened!\n";
		std::cout << "Possible errors:\n";
		std::cout << "1. The file does not exist.\n";
		std::cout << "2. The path was not found.\n";
		std::exit(-1);
	}
	setfile.close();
}

// @brief help message to run the software
auto HelpMessage(char* a) -> void {
	std::cerr << "Pay attention to the command line arguments!\n"
				<< "Usage: " << a << " [-h] [-g GENERAL OPTIONS] [-s SPECIES] [-t TOY]"
        << "\nOptions:\n"
        << "\t-h,--help\n\tShow this help message\n"
        << "\t-g,--general-options\n\tSpecify the file containing PDFs, data, output rootfiles paths and other info.\n"
        << "\t-s,--species-list\n\tSpecify the file containing the number of parameters, the list of species (also if they are free/fixed/constrained) and the min/man energy range\n"
		<< "\t-o,--output\n\tSpecify the path containing the output, without the extension since the fitter will create an output.root and output.txt.\n"
        << "\t-t,--toy-rates\n\tTO BE WRITTEN"
		<< std::endl;
}

namespace CMDLParser {

auto Parse(int argc, char* argv[]) -> NuFitCmdlArgs{

	auto args = std::make_unique<NuFitCmdlArgs>();

	if(argc < 7 || (std::strcmp(argv[1], "--help")) == 0 || (std::strcmp(argv[1], "-h") == 0)){
		NuFitter::HelpMessage(argv[0]);
		std::exit(-1);
	}
	for (int i = 1; i < argc-1; i++){
			if (std::strcmp(argv[i], "--general-options") == 0 || std::strcmp(argv[i], "-g") == 0){
				args->gen = argv[i+1];
				i++;
			}
			else if (std::strcmp(argv[i], "--species-list") == 0 || std::strcmp(argv[i], "-s") == 0){
				args->spec = argv[i+1];
				i++;
			}
			else if (std::strcmp(argv[i], "--toy-rates") == 0 || std::strcmp(argv[i], "-t") == 0){
				args->toy = argv[i+1];
				i++;
			}
			else if (std::strcmp(argv[i], "--output") == 0 || std::strcmp(argv[i], "-o") == 0){
				args->output_name = argv[i+1];
				i++;
			}
			else {
				std::cout << "ERROR: something went wrong. Please read carefully the instructions below.\n";
				NuFitter::HelpMessage(argv[0]);
				std::exit(-1);
		}
	}

	return *args;
}
}  // namespace CMDLParser


namespace ConfigParser {

auto ParseGenOpts(std::unique_ptr<NuFitConfig> &config, std::string filename) -> void {
	// Make sure the file is valid
	NuFitter::ErrorReading(filename);

	// Loop through the file
	std::string line, key, value;
	std::ifstream ReadGen;
	ReadGen.open(filename);

	while (std::getline(ReadGen, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // Read from line
		std::stringstream line_stream(line);

		// Ignore commented lines
		line_stream >> key;
		if (std::strncmp(key.c_str(), "#", 1) == 0) {
			continue;
		}

		// Continue reading line
		line_stream >> value;

		// Assign the values to their keys
		if (key.find("Hist") != std::string::npos) {
			config->data_hist_names.push_back(value);
		}
		else if (key == "PDFsRootfile") config->pdf_name = value;
		else if (key == "DataRootfile") config->data_name = value;
		else if (key == "Lifetime") config->lifetime = std::stod(value);
		else if (key == "MassTarget") config->mass_target = std::stod(value);
		else if (key == "emax") config->emax = std::stod(value);
		else if (key == "emin") config->emin = std::stod(value);
		else if (key == "ToyData") config->doToyData_ = std::stoi(value);
		else if (key == "Hesse") config->doHesse = std::stoi(value);
		else if (key == "Minos") config->doMinos = std::stoi(value);
		else if (key == "Likelihood") config->likelihood = value;
		// Todo error handling in case some values are missing
	}
	config->exposure = config->lifetime * config->mass_target;

    std::cout << "Found " << config->data_hist_names.size()
              << " histograms:" << std::endl;
    for (auto i = 1U; i <= config->data_hist_names.size(); i++) {
        std::cout << "\t" << i << ": " << config->data_hist_names[i-1] << std::endl;
    }
}

// @brief Parse the species_list content and save to config
auto ParseSpeciesList(std::unique_ptr<NuFitConfig>& config,
	                  std::string filename) -> void {
	// Make sure the file is valid
	NuFitter::ErrorReading(filename);

	// Initialise variables
	std::string line, word;
	auto nPDFs{0}, nParams{0};

	// Open the file
	std::ifstream ReadSpec;
	ReadSpec.open(filename);

	// For each line in the file
	while (std::getline(ReadSpec, line)) {
		std::stringstream line_stream(line);
		auto nElements {0U};

		while (line_stream >> word) {
			// Ignore commented lines
			if (nElements == 0 && std::strncmp(word.c_str(), "#", 1) == 0) {
				break;
			}

			// Read the (expected) variables
			switch (nElements) {
				case 0:
					nPDFs += 1;
					config->pdf_names.push_back(word);
					break;
				case 1:
					// Check if param_names already contains the same word
					if (std::find(config->param_names.begin(), config->param_names.end(), word)
						== config->param_names.end()) {
						nParams++;
					}
					config->param_names.push_back(word);
					break;
				case 2:
					config->param_initial_guess.push_back(std::stod(word)*config->exposure);
					break;
				case 3:
					config->param_lowerlim.push_back(std::stod(word)*config->exposure);
					break;
				case 4:
					config->param_upperlim.push_back(std::stod(word)*config->exposure);
					break;
				case 5:
					config->param_stepsize.push_back(std::stod(word));
					break;
				case 6:
					config->param_fixed.push_back(std::stoi(word));
					// Fixed parameters do not count in the fit
					if (std::stoi(word) == 1) nParams--;
					break;
				case 7:
					config->hist_id.push_back(std::stoi(word));
					break;
				case 8:
					config->param_eff.push_back(std::stod(word));
					break;
				default:
					std::cout << "[Warning] In file: " + filename +
						": Too many arguments in a line." << std::endl;
			}
			nElements++;
		}
		// Make sure line was correctly parsed
		if (nElements == 0) continue;
		if (nElements < 8) {
			throw std::invalid_argument("A line in file: '" + filename +
				"' does not have enough arguments to be valid.."
				+ std::to_string(nElements));
		}
	}
	ReadSpec.close();
    // TODO: additional error handling

	// Count how many species exist in each histogram
	auto max = *std::max_element(std::begin(config->hist_id),
	                             std::end(config->hist_id));
	std::vector<unsigned int> tmp_nSp(max);
	for (auto el : config->hist_id) {
		tmp_nSp[el-1]++;
	}
	config->nSp_histos = tmp_nSp;

	// Bookkeeping
	config->npdfs = nPDFs;
	config->nparams = nParams;

}

auto Parse(NuFitCmdlArgs args) -> NuFitConfig {
	auto config = std::make_unique<NuFitConfig>();
	config->output_name = args.output_name;

	// -------------------------------------------------------------------------
	// Read the general_options config file
	ParseGenOpts(config, args.gen);

	// -------------------------------------------------------------------------
	// Read the species-list
	ParseSpeciesList(config, args.spec);

	// -------------------------------------------------------------------------
	// TODO: make sure pdfs are compatible
	// Get nbins from the data hists or from the PDFs if we want to use ToyData
	if(config->doToyData_ == false){
		for (auto i = 0U; i < config->data_hist_names.size(); i++) {
			// Only load the histogram if it is used!
			if (std::find(config->hist_id.begin(), config->hist_id.end(), i+1)
		    == config->hist_id.end()) continue;
				TFile *fdata = new TFile(config->data_name.c_str());
				TH1D* hdata = (TH1D*)fdata->Get(config->data_hist_names[i].c_str());
				config->nbins.push_back(hdata->GetNbinsX());
				fdata->Close();
	}
} else {
	for (auto i = 0U; i < config->data_hist_names.size(); i++) {
		TFile *fpdf = new TFile(config->pdf_name.c_str());
		TH1D* hpdf = (TH1D*)fpdf->Get(config->pdf_names[0]);
		config->nbins.push_back(hpdf->GetNbinsX());
		fpdf->Close();
	}
}


	return *config;
}

}  // namespace ConfigParser
}  // namespace NuFitter
