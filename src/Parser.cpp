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
	
	if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1],"-h") == 0){
		std::cerr << "Usage: " << argv[0] << " [-h] [-g GENERAL OPTIONS] [-s SPECIES] [-t TOY]"
              << "Options:\n"
              << "\t-h,--help\tShow this help message\n"
              << "\t-g,--gen GENERAL OPTIONS\tSpecify the file containing PDFs, data, output rootfiles paths and other info.\n"
              << "\t-s,--species SPECIES\tSpecify the file containing the number of parameters, the list of species (also if they are free/fixed/constrained) and the min/man energy range\n"
              << "\t-t,--toy TOY\tTO BE WRITTEN"
              << std::endl;

		exit(-1);
	}

	for (int i = 1; i < argc; i++){
		if (i + 1 != argc){
			if (strcmp(argv[i], "--gen") == 0 || strcmp(argv[i], "-g") == 0){
				args->gen = argv[i+1];
				i++;
			} else 
				if (strcmp(argv[i], "--species") == 0 || strcmp(argv[i], "-s") == 0){
                        	        args->spec = argv[i+1];
                                	i++;
                        	}
			 else 
				if (strcmp(argv[i], "--toy") == 0 || strcmp(argv[i], "-t") == 0){
                                	args->toy = argv[i+1];
                               		i++;
				}
			 else {	
				 std::cout << "ERROR: something went wrong. To find out more, please try: " << argv[0] << " -h.\n";
				 exit(-1);
			 }

		}
	}

	return *args;
}
}  // namespace CMDLParser

namespace ConfigParser {

auto Parse(NuFitCmdlArgs args) -> NuFitConfig {

    // Set-up ifstream from config file
	std::ifstream ReadGen, ReadSpec;
	ReadGen.open(args.gen.c_str());
	ReadSpec.open(args.spec.c_str());

	if(ReadGen.fail()){
		std::cout << "Opening " << args.gen.c_str() << " for reading.\n";
		std::cout <<"The "<< args.gen.c_str() <<" file could not be opened!\n";
		std::cout << "Possible errors:\n";
		std::cout <<"1. The file does not exist.\n";
		std::cout <<"2. The path was not found.\n";
		exit(-1);
	}

	if(ReadSpec.fail()){
                std::cout << "Opening " << args.spec.c_str() << " for reading.\n";
                std::cout <<"The "<< args.spec.c_str() <<" file could not be opened!\n";
                std::cout << "Possible errors:\n";
                std::cout <<"1. The file does not exist.\n";
                std::cout <<"2. The path was not found.\n";
                exit(-1);
        }

	auto config = std::make_unique<NuFitConfig>();

	//quite rough, to be changed
    // Define placeholder variables to store the values
	std::string output, datafile, pdffile, datahist, toy, hesse, minos;
	double min, max;
	std::string appo1;
	double appo2;
	int npar = 0;

    // Read the file
    ReadGen >> output >> pdffile >> datafile >> datahist >> toy >> hesse >> minos >> min >> max;
 	
    // Save and write to new NuFitConfig object
	std::string unused;
	while(std::getline(ReadSpec,unused))
   	++npar;

	std::cout << std::endl << "npar " << npar << std::endl;
	ReadSpec.close();
	ReadSpec.open(args.spec.c_str());

	for(int i = 1; i <= npar; i++){
		
		ReadSpec >> appo1;
		config->param_names.push_back(appo1);	
		ReadSpec >> appo2;
		config->param_initial_guess.push_back(appo2);
		ReadSpec >> appo2;
		config->param_lowerlim.push_back(appo2);
		ReadSpec >> appo2;
		config->param_upperlim.push_back(appo2);
		ReadSpec >> appo2;
		config->param_stepsize.push_back(appo2);
	}	

	config->nparams = npar;
	config->npdfs = npar;
	config->output_name = output;
    	config->data_name = datafile;
    	config->pdf_name = pdffile;
   	config->histo_data = datahist;
	config->emin = min;
	config->emax = max;

	if(toy == "no")	config->doToyData_ = false;
	else	config->doToyData_ = true;

	if(hesse == "no") config->doHesse = false;
        else    config->doHesse = true;

	if(minos == "no") config->doMinos = false;
        else    config->doMinos = true;



	// Get nbins from the data hist
	// TODO: make sure pdfs are compatible?
	TFile *fdata = new TFile(config->data_name.c_str());
	TH1D* hdata = (TH1D*)fdata->Get(config->histo_data.c_str());
	config->nbins = hdata->GetNbinsX();
	fdata->Close();

	ReadGen.close();
	ReadSpec.close();

	return *config;
}

}  // namespace ConfigParser
}  // namespace NuFitter
