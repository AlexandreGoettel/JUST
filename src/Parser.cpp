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
#include <cstring>
#include <cstdlib>
// ROOT includes
#include "TFile.h"
#include "TH1D.h"
// Project includes
#include "Parser.h"

//============================================================================
// Implementations

namespace NuFitter{

namespace CMDLParser {

auto Parse(int argc, char* argv[]) -> NuFitCmdlArgs{

	auto args = std::make_unique<NuFitCmdlArgs>();

	if(argc < 4 || (std::strcmp(argv[1], "--help")) == 0 || (std::strcmp(argv[1], "-h") == 0)){
		NuFitter::HelpMessage(argv[0]);
		std::exit(-1);
	}
	for (int i = 1; i < argc; i++){
		if (i + 1 != argc){
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
	}

	return *args;
}
}  // namespace CMDLParser


namespace ConfigParser {

auto Parse(NuFitCmdlArgs args) -> NuFitConfig {
	auto config = std::make_unique<NuFitConfig>();
	config->output_name = args.output_name;

	// -------------------------------------------------------------------------
	// Read the general_options config file
	NuFitter::ErrorReading(args.gen);

	std::vector<std::string> labels = {"PDFsRootfile", "DataRootfile",
	    "HistoName", "Lifetime", "MassTarget", "emin", "emax", "ToyData",
		"Hesse", "Minos", "Likelihood"};

	config->pdf_name = NuFitter::GetValue(labels.at(0),args.gen);
	config->data_name = NuFitter::GetValue(labels.at(1),args.gen);
	config->histo_data = NuFitter::GetValue(labels.at(2),args.gen);
	config->lifetime = std::stod(NuFitter::GetValue(labels.at(3),args.gen));
	config->mass_target = std::stod(NuFitter::GetValue(labels.at(4),args.gen));
	config->emin = std::stod(NuFitter::GetValue(labels.at(5),args.gen));
	config->emax = std::stod(NuFitter::GetValue(labels.at(6),args.gen));
	std::istringstream(NuFitter::GetValue(labels.at(7),args.gen)) >> config->doToyData_;
	std::istringstream(NuFitter::GetValue(labels.at(8),args.gen)) >> config->doHesse;
	std::istringstream(NuFitter::GetValue(labels.at(9),args.gen)) >> config->doMinos;
	config->likelihood = NuFitter::GetValue(labels.at(10),args.gen);

	config->exposure = config->lifetime * config->mass_target;

	// -------------------------------------------------------------------------
	// Read the species-list
	std::ifstream ReadSpec;
	ReadSpec.open(args.spec);

	std::string line, word;
	auto nSpecies{0};
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
					config->param_names.push_back(word);
					break;
				case 1:
					config->param_initial_guess.push_back(std::stod(word)*config->exposure);
					break;
				case 2:
					config->param_lowerlim.push_back(std::stod(word)*config->exposure);
					break;
				case 3:
					config->param_upperlim.push_back(std::stod(word)*config->exposure);
					break;
				case 4:
					config->param_stepsize.push_back(std::stod(word));
					break;
				case 5:
					config->param_fixed.push_back(std::stoi(word));
					break;
				case 6:
					config->param_constr_sigma.push_back(std::stod(word));
					break;
				default:
					std::cout << "[Warning] In file: " + args.spec +
						": Too many arguments in a line." << std::endl;
			}
			nElements++;
		}
		// Make sure line is correct
		if (nElements == 0) continue;
		if (nElements < 6) {
			throw std::invalid_argument("A line in file: '" + args.spec +
				"' does not have enough arguments to be valid.."
				+ std::to_string(nElements));
		}
		if (nElements == 6) {  // If no param constraint was given, fill zero
			config->param_constr_sigma.push_back(0.);
		}

		// Bookkeeping
		nSpecies++;
	}

	config->nparams = nSpecies;
	config->npdfs = nSpecies;
	ReadSpec.close();
	// -------------------------------------------------------------------------

	// Get nbins from the data hist
	// TODO: make sure pdfs are compatible?
	TFile *fdata = new TFile(config->data_name.c_str());
	TH1D* hdata = (TH1D*)fdata->Get(config->histo_data.c_str());
	config->nbins = hdata->GetNbinsX();
	fdata->Close();

	return *config;
}

}  // namespace ConfigParser

// @brief problems in opening or reading input files
auto ErrorReading(std::string& filename) -> void {
	std::ifstream setfile(filename);
	if(setfile.fail()){
		std::cout << "Opening " << filename << " for reading.\n";
		std::cout <<"The "<< filename <<" file could not be opened!\n";
		std::cout << "Possible errors:\n";
		std::cout <<"1. The file does not exist.\n";
		std::cout <<"2. The path was not found.\n";
		std::exit(-1);
	}
	setfile.close();
}

// @brief help message to run the software
auto HelpMessage(char* a) -> void {
	std::cerr << "Usage: " << a << " [-h] [-g GENERAL OPTIONS] [-s SPECIES] [-t TOY]"
        << "\nOptions:\n"
        << "\t-h,--help\n\tShow this help message\n"
        << "\t-g,--general-options\n\tSpecify the file containing PDFs, data, output rootfiles paths and other info.\n"
        << "\t-s,--species-list\n\tSpecify the file containing the number of parameters, the list of species (also if they are free/fixed/constrained) and the min/man energy range\n"
		<< "\t-o,--output\n\tSpecify the path containing the output, without the extension since the fitter will create an output.root and output.txt.\n"
        << "\t-t,--toy-rates\n\tTO BE WRITTEN"
		<< std::endl;
}

// @brief Read and Fill for general_options.txt
template<class T>
auto ReadAndFill_Gen(std::ifstream& filename, T& var1, T& var2) -> void {
	std::string appo;
	filename >> appo;  // read the labels
	filename >> var1;
	var2 = var1;
}
// @brief Loop in filename.txt, search for "variable" and return its value
inline auto GetValue(std::string& variable, std::string& filename) -> std::string {
	std::ifstream setfile(filename);
	std::string val;
	std::string var;
	bool anyfound(false);

	while(1){
		setfile >> var >> val;
		if (!setfile.good()) break;
		if(var == variable){
			anyfound = true;
			break;
		}
	}

	if(!anyfound) std::cout << "warning!!: I didn't find " << variable
	                        << "...setting it to 0" << std::endl;

	setfile.close();
	return val;
}

}  // namespace NuFitter
