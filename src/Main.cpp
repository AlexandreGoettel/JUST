//! @file	src/Main.cpp
//! @brief	Declaration of main guideline functions
//! @author	Alexandre Göttel
//! @date	2021-03-18

//============================================================================
// Standard includes
// #include <memory>
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
	const NuFitConfig *config = ConfigParser::Parse(cmdl_args);

	// Read the PDFs to fit to the data
	NuFitPDFs *pdfs = PDFs::Read(config);

	// Perform the fit with real or toy data
	NuFitResults *results;
	if (config->doToyData_) {
		// Generate toy data for the fit
		NuFitToyData *data = generateToyData(config);
		results = MCFit::Fit(data, pdfs, config);
	} else {
		// Read the data histogram(s)
		NuFitData *data = Data::Read(config);
		results = MCFit::Fit(data, pdfs, config);
	}

	// Save/plot results
	ProcessResults(results);
}
