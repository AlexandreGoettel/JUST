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

	if(argc < 4 || (strcmp(argv[1], "--help")) == 0 || (strcmp(argv[1], "-h") == 0)){
		NuFitter::HelpMessage(argv[0]);
		exit(-1);
	}
	for (int i = 1; i < argc; i++){
		if (i + 1 != argc){
			if (strcmp(argv[i], "--gen") == 0 || strcmp(argv[i], "-g") == 0){
				args->gen = argv[i+1];
				i++;
			}
			else 
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
						 std::cout << "ERROR: something went wrong. Please read carefully the instructions below.\n";
						 NuFitter::HelpMessage(argv[0]);
						 exit(-1);
		 			}
			}
	}

		

	return *args;
}
}  // namespace CMDLParser


namespace ConfigParser {

auto Parse(NuFitCmdlArgs args) -> NuFitConfig {

	auto config = std::make_unique<NuFitConfig>();

	std::ifstream ReadGen, ReadSpec;

	ReadGen.open(args.gen);
	NuFitter::ErrorReading(ReadGen,args.gen);

	//placeholder variables
	std::string output, pdffile, datafile, datahist;
        bool toy, hesse, minos = false;
        double ltime, mass, min, max = 0;

	//general_options.txt: read and fill the NuFitConfig variable
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

	ReadGen.close();

	//count the number of species from species_list.txt
	int N = HowManySpecies(ReadSpec,args.spec);
	config->nparams = N;
	config->npdfs = N;

	ReadSpec.open(args.spec);
	NuFitter::ErrorReading(ReadSpec,args.spec);

	//placeholder variables
	TString namepar;
	double inguess, lowlim, uplim, step = 0;
	
	//loop over the number of species and fill the NuFitConfig variable
	for(int i = 0; i < N; i++){

		NuFitter::ReadAndFill_Spec(ReadSpec,namepar,config->param_names);
		NuFitter::ReadAndFill_Spec(ReadSpec,inguess,config->param_initial_guess);
		NuFitter::ReadAndFill_Spec(ReadSpec,lowlim,config->param_lowerlim);
		NuFitter::ReadAndFill_Spec(ReadSpec,uplim,config->param_upperlim);
		NuFitter::ReadAndFill_Spec(ReadSpec,step,config->param_stepsize);
		
	}
	
	ReadSpec.close();

	// Get nbins from the
	//
	// data hist
	// TODO: make sure pdfs are compatible?
	TFile *fdata = new TFile(config->data_name.c_str());
	TH1D* hdata = (TH1D*)fdata->Get(config->histo_data.c_str());
	config->nbins = hdata->GetNbinsX();
	fdata->Close();

	return *config;
}

}  // namespace ConfigParser

//problems in opening or reading input files
auto ErrorReading(const std::ifstream& filename, const std::string& s) -> void {

        if(filename.fail()){
                std::cout << "Opening " << s << " for reading.\n";
                std::cout <<"The "<< s <<" file could not be opened!\n";
                std::cout << "Possible errors:\n";	   
		std::cout <<"1. The file does not exist.\n";
                std::cout <<"2. The path was not found.\n";
                exit(-1);
        }
}

//help message to run the software
auto HelpMessage(char* a) -> void {

	std::cerr << "Usage: " << a << " [-h] [-g GENERAL OPTIONS] [-s SPECIES] [-t TOY]"
        << "\nOptions:\n"
        << "\t-h,--help\tShow this help message\n"
        << "\t-g,--gen GENERAL OPTIONS\tSpecify the file containing PDFs, data, output rootfiles paths and other info.\n"
        << "\t-s,--species SPECIES\tSpecify the file containing the number of parameters, the list of species (also if they are free/fixed/constrained) and the min/man energy range\n"
        << "\t-t,--toy TOY\tTO BE WRITTEN"
	<< std::endl;
}

//Read and Fill for general_options.txt
template<class T> auto ReadAndFill_Gen(std::ifstream& filename, T& var1, T& var2) -> void {

	std::string appo;
	filename >> appo;//read the labels
	filename >> var1;
	var2 = var1;

}

//Read and Fill for species_list.txt
template<class T> auto ReadAndFill_Spec(std::ifstream& filename, T& var1, std::vector<T>& var2) -> void {

        filename >> var1;
        var2.push_back(var1);

}

//count the number of species from species_list.txt
auto HowManySpecies(std::ifstream& filename, const std::string& s) -> int {

	filename.open(s);
	std::string unused;
	int n = 0;

        while(std::getline(filename,unused))	++n;
	
	filename.close();

	return n;

}
}  // namespace NuFitter
