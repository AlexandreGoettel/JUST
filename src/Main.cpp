//! @file	src/Main.cpp
//! @brief	Declaration of main guideline functions
//! @author	Alexandre Göttel
//! @date	2021-03-18

//============================================================================
// Standard includes
#include <memory>
#include <vector>
#include <iostream>
// Project includes
#include "Parser.h"
#include "DataReader.h"
#include "ToyDataGenerator.h"  // Insert if here?
#include "Fitter.h"
#include "OutputManager.h"

//============================================================================
// Method implementation

//____________________________________________________________________________
//! @param  argc  Command line (CMDL) argument count.
//! @param  argv  Command line (CMDL) argument values.
//! @return Status code, which is 0 if there were no problems.
NuFitConfig const *config;
auto main(int argc, char* argv[]) -> int {
	using namespace NuFitter;

	// Parse CMDL arguments
	const NuFitCmdlArgs cmdl_args = CMDLParser::Parse(argc, argv);

	// Parse config files
	// NuFitConfig const *config = new NuFitConfig(cmdl_args);
	config = new NuFitConfig(cmdl_args);

	// Read the PDFs to fit to the data
	NuFitPDFs *pdfs = new NuFitPDFs();
	pdfs->Read(config->pdf_names, config->npdfs, config->hist_id);

	// Perform the fit with real or toy data
	if (config->ToyData != 0) {
		// Read the pdfs used to generate the toy data
		NuFitPDFs *pdfs_toy = new NuFitPDFs();
		pdfs_toy->Read(config->pdf_names_toy, config->npdfs_toy, config->hist_id_toy);

		// Generate toy data for the fit
		NuFitToyData *data = new NuFitToyData(pdfs_toy);

		// Fit the toy-data
		auto results = MCFit::Fit(data, pdfs);
		ProcessResults(data, pdfs_toy, pdfs, results);

		delete data;
		delete pdfs_toy;
	} else {
		// Read the data histogram(s)
		NuFitData *data = new NuFitData();
		data->Read();

		// Fit the data
		auto results = MCFit::Fit(data, pdfs);
		ProcessResults(data, pdfs, results);

		delete data;
	}
	delete pdfs;
	delete config;
}
