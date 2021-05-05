//! @file	src/Parser.cpp
//! @brief	Implementation of classes to parse cmdl args and config files
//! @author	Alexandre Göttel
//! @date	2021-03-18

//============================================================================
// Includes
// Standard includes
#include <memory>
#include <fstream>
#include <iostream>
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

	// -------------------------------------------------------------------------
	// Read the general_options config file
	std::ifstream ReadGen;
	ReadGen.open(args.gen);
	NuFitter::ErrorReading(ReadGen,args.gen);

	//placeholder variables
	std::string output, pdffile, datafile, datahist, likelihood;
	bool toy, hesse, minos = false;
	double ltime, mass, min, max = 0;

	// general_options.txt: read and fill the NuFitConfig variable
	// TODO: loop over the file, can probably auto placeholder
	NuFitter::ReadAndFill_Gen(ReadGen,output,config->output_name);
	NuFitter::ReadAndFill_Gen(ReadGen,pdffile,config->pdf_name);
	NuFitter::ReadAndFill_Gen(ReadGen,datafile,config->data_name);
	NuFitter::ReadAndFill_Gen(ReadGen,datahist,config->histo_data);
	NuFitter::ReadAndFill_Gen(ReadGen,ltime,config->lifetime);
	NuFitter::ReadAndFill_Gen(ReadGen,mass,config->mass_target);
	NuFitter::ReadAndFill_Gen(ReadGen,min,config->emin);
	NuFitter::ReadAndFill_Gen(ReadGen,max,config->emax);
	NuFitter::ReadAndFill_Gen(ReadGen,toy,config->doToyData_);
	NuFitter::ReadAndFill_Gen(ReadGen,hesse,config->doHesse);
	NuFitter::ReadAndFill_Gen(ReadGen,minos,config->doMinos);
	NuFitter::ReadAndFill_Gen(ReadGen,likelihood,config->likelihood);
	ReadGen.close();


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
					config->param_initial_guess.push_back(std::stod(word));
					break;
				case 2:
					config->param_lowerlim.push_back(std::stod(word));
					break;
				case 3:
					config->param_upperlim.push_back(std::stod(word));
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
auto ErrorReading(const std::ifstream& filename, const std::string& s) -> void {
	if(filename.fail()){
		std::cout << "Opening " << s << " for reading.\n";
		std::cout <<"The "<< s <<" file could not be opened!\n";
		std::cout << "Possible errors:\n";
		std::cout <<"1. The file does not exist.\n";
		std::cout <<"2. The path was not found.\n";
		std::exit(-1);
	}
}

// @brief help message to run the software
auto HelpMessage(char* a) -> void {
	std::cerr << "Usage: " << a << " [-h] [-g GENERAL OPTIONS] [-s SPECIES] [-t TOY]"
        << "\nOptions:\n"
        << "\t-h,--help\n\tShow this help message\n"
        << "\t-g,--general-options\n\tSpecify the file containing PDFs, data, output rootfiles paths and other info.\n"
        << "\t-s,--species-list\n\tSpecify the file containing the number of parameters, the list of species (also if they are free/fixed/constrained) and the min/man energy range\n"
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

}  // namespace NuFitter
