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
#include "TMath.h"
#include "TRandom3.h"
// Project includes
#include "Parser.h"

#define PI 3.141592653589793

//============================================================================
// Implementations

namespace NuFitter {

// @brief Get the index of element el in vector v (first occurence)
// TODO: move to more sensible place to share with fitter?
template <class T>
auto getIndexOf_(T el, std::vector<T> v) -> int {
	auto it = std::find(v.begin(), v.end(), el);

	// If element was found, calculate the index
	if (it != v.end()) {
		return it - v.begin();
	}
	else {
		return -1;
	}
}

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
}  // namespace NuFitter

auto NuFitConfig::ParseGenOpts(std::string filename) -> void {
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
			data_hist_names.push_back(value);
		}
		else if (key == "PDFsRootfile") pdf_name = value;
		else if (key == "DataRootfile") data_name = value;
		else if (key == "Lifetime") lifetime = std::stod(value);
		else if (key == "LSDensity_g/mL") densityLS = std::stod(value);
		else if (key == "Radius_m") radius = std::stod(value);
		else if (key == "DAQTime") daq_time = std::stod(value);
		else if (key == "Efficiency") efficiency = std::stod(value);
		else if (key == "emax") emax = std::stod(value);
		else if (key == "emin") emin = std::stod(value);
		else if (key == "ToyData") ToyData = std::stoi(value);
		else if (key == "Hesse") doHesse = std::stoi(value);
		else if (key == "Minos") doMinos = std::stoi(value);
		else if (key == "Likelihood") likelihood = value;
		else if (key == "TargetMass") mass_target = std::stod(value);
		else if (key == "Seed") seed = std::stoi(value);
		// Todo error handling in case some values are missing
	}

	// Give the option to define the target mass directly, or by geometry
	if (mass_target <= 0.) {
		std::cout << "calculating target mass from geometry..." << std::endl;
		mass_target = 4 * PI * pow(radius,3) * densityLS / 3000.;
	}
	exposure = lifetime * mass_target * daq_time * efficiency;
	assert(exposure > 0);

    std::cout << "Found " << data_hist_names.size()
              << " histograms:" << std::endl;
    for (auto i = 1U; i <= data_hist_names.size(); i++) {
        std::cout << "\t" << i << ": " << data_hist_names[i-1] << std::endl;
    }
}

// @brief Parse the species_list content and save to config
auto NuFitConfig::ParseSpeciesList(std::string filename) -> void {
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
					pdf_names.push_back(word);
					break;
				case 1:
					// Check if param_names already contains the same word
					if (std::find(param_names.begin(), param_names.end(), word)
						== param_names.end()) {
						nParams++;
					}
					param_names.push_back(word);
					break;
				case 2:
					param_initial_guess.push_back(std::stod(word)*exposure);
					break;
				case 3:
					param_lowerlim.push_back(std::stod(word)*exposure);
					break;
				case 4:
					param_upperlim.push_back(std::stod(word)*exposure);
					break;
				case 5:
					param_stepsize.push_back(std::stod(word)*exposure);
					break;
				case 6:
					param_fixed.push_back(std::stoi(word));
					// Fixed parameters do not count in the fit
					if (std::stoi(word) == 1) nParams--;
					break;
				case 7:
					hist_id.push_back(std::stoi(word));
					break;
				case 8:
					param_eff.push_back(std::stod(word));
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
	// TODO: count the number of different hists instead. Could be that
	// people use hist 1 and 3 but not 2
	auto max = *std::max_element(std::begin(hist_id),
	                             std::end(hist_id));
	std::vector<unsigned int> tmp_nSp(max);
	for (auto el : hist_id) {
		tmp_nSp[el-1]++;
	}
	nSp_histos = tmp_nSp;

	// Bookkeeping
	npdfs = nPDFs;
	nparams = nParams;
}

// @brief Parse the toy_rates content and save to config
auto NuFitConfig::ParseToyRates(std::string filename) -> void {
	// Make sure the file is valid
	NuFitter::ErrorReading(filename);

	// Initialise variables
	std::string line, word;
	auto nPDFs{0};

	// Open the file
	std::ifstream ReadToy;
	ReadToy.open(filename);

	// For each line in the file
	while (std::getline(ReadToy, line)) {
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
					pdf_names_toy.push_back(word);
					break;
				case 1:
					param_initial_guess_toy.push_back(std::stod(word)*exposure);
					break;
				case 2:
					hist_id_toy.push_back(std::stoi(word));
					break;
				case 3:
					param_eff_toy.push_back(std::stod(word));
					break;
				default:
					std::cout << "[Warning] In file: " + filename +
						": Too many arguments in a line." << std::endl;
			}
			nElements++;
		}
		// Make sure line was correctly parsed
		if (nElements == 0) continue;
		if (nElements < 3) {
			throw std::invalid_argument("A line in file: '" + filename +
				"' does not have enough arguments to be valid.."
				+ std::to_string(nElements));
		}
	}
	ReadToy.close();
    // TODO: additional error handling

	// Count how many species exist in each histogram
	auto max = *std::max_element(std::begin(hist_id_toy),
	                             std::end(hist_id_toy));
	std::vector<unsigned int> tmp_nSp(max);
	for (auto el : hist_id_toy) {
		tmp_nSp[el-1]++;
	}
	nSp_histos_toy = tmp_nSp;

	// Bookkeeping
	npdfs_toy = nPDFs;
}

// @brief Constructor for NuFitConfig
NuFitConfig::NuFitConfig(NuFitter::NuFitCmdlArgs args) {
	output_name = args.output_name;
	// -------------------------------------------------------------------------
	// Read the general_options config file
	ParseGenOpts(args.gen);

	// -------------------------------------------------------------------------
	// Read the species-list
	ParseSpeciesList(args.spec);

	// -------------------------------------------------------------------------
	// Read the toy-rates
	if(ToyData != 0){
		ParseToyRates(args.toy);

		std::cout << "Running random number generator for toy-data..." << std::endl;
		// Sample species counts here to pass to toyDataGenerator
		gRandom = new TRandom3(0);
		gRandom->SetSeed(seed);
		for (auto t = 0U; t < ToyData; t++) {  // For each toy data fit
			std::vector<unsigned long int> current_sampled_counts;
			for (auto i = 0U; i < npdfs_toy; i++) {  // For each pdf
				auto n_expected = param_initial_guess_toy[i]*param_eff_toy[i];
				unsigned long n_sampled;
				// Use Gaus() for large numbers to avoid numeric int limits
				if (n_expected > 10000) {
					n_sampled = gRandom->Gaus(n_expected, sqrt(n_expected));
				} else {
					n_sampled = gRandom->Poisson(n_expected);
				}
				current_sampled_counts.push_back(n_sampled);
			}
			param_sampled.push_back(current_sampled_counts);
		}

		// Create a paramVector object for the toy data parameters
		std::vector<TString> used_names, used_names_fixed;
		for (auto i = 0U; i < npdfs_toy; i++) {
			NuFitter::paramData current_paramData {i, hist_id_toy[i]};
			auto name = pdf_names_toy[i];

			if (std::find(used_names.begin(), used_names.end(), name) == used_names.end()) {
				std::vector<NuFitter::paramData> tmp_paramVector;
				tmp_paramVector.push_back(current_paramData);
				paramVector_toy.push_back(tmp_paramVector);
				used_names.push_back(name);
			} else {
				auto idx_name = NuFitter::getIndexOf_(name, used_names);
				assert(idx_name != -1);
				paramVector_toy[idx_name].push_back(current_paramData);
			}
		}
	}

	// -------------------------------------------------------------------------
	// TODO: make sure pdfs are compatible
	// Get nbins from the data hists or from the PDFs if we want to use ToyData
	if(ToyData == 0){
		for (auto i = 0U; i < data_hist_names.size(); i++) {
			// Only load the histogram if it is used!
			if (std::find(hist_id.begin(), hist_id.end(), i+1)
			    == hist_id.end()) continue;
			TFile *fdata = new TFile(data_name.c_str());
			TH1D* hdata = (TH1D*)fdata->Get(data_hist_names[i].c_str());
			nbins.push_back(hdata->GetNbinsX());
			fdata->Close();
		}
	} else {
		for (auto i = 0U; i < data_hist_names.size(); i++) {
			// Only load the histogram if it is used!
			if (std::find(hist_id_toy.begin(), hist_id_toy.end(), i+1)
			    == hist_id_toy.end()) continue;
			TFile *fpdf = new TFile(pdf_name.c_str());
			TH1D* hpdf = (TH1D*)fpdf->Get(pdf_names[0]);
			nbins.push_back(hpdf->GetNbinsX());
			fpdf->Close();
		}
	}
}

// @brief Destructor for NuFitConfig
NuFitConfig::~NuFitConfig() {}
