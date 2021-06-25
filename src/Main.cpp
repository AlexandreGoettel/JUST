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
auto main(int argc, char* argv[]) -> int {
	using namespace NuFitter;

	// Parse CMDL arguments
	const NuFitCmdlArgs cmdl_args = CMDLParser::Parse(argc, argv);

	// Parse config files
	const NuFitConfig config = ConfigParser::Parse(cmdl_args);

	// Read the PDFs to fit to the data
	NuFitPDFs *pdfs = PDFs::Read(config);

	// Perform the fit with real or toy data
	if (config.ToyData != 0) {
		// Read the pdfs used to generate the toy data
		NuFitPDFs *pdfs_toy = Toy::Read(config);
		// Generate toy data for the fit
		auto data = ToyData::Initialise(config, pdfs_toy);
		auto results = MCFit::Fit(data, pdfs, config);
		ProcessResults(data, pdfs, config, results);
	} else {
		// Read the data histogram(s)
		auto data = Data::Read(config);
		auto results = MCFit::Fit(data, pdfs, config);
		ProcessResults(data, pdfs, config, results);
	}
}
