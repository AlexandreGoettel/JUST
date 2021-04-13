//! @file	include/Fitter.h
//! @brief	Declaration of classes to parse cmdl args and config files
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef PARSER_H_
#define PARSER_H_

//============================================================================
// Includes
#include <vector>
#include "TString.h"

//============================================================================
// Class definition

namespace NuFitter {
class NuFitCmdlArgs {
public:
	NuFitCmdlArgs() = default;  // constructor
	~NuFitCmdlArgs() = default;  // destructor
	NuFitCmdlArgs(const NuFitCmdlArgs&) = default;  // copy constructor
	NuFitCmdlArgs(NuFitCmdlArgs&&) = default;  // move constructor
	NuFitCmdlArgs &operator=(const NuFitCmdlArgs&) = default;  // copy assignment
	NuFitCmdlArgs &operator=(NuFitCmdlArgs&&) = default;  // move assignment
};

class NuFitConfig {
public:
	NuFitConfig();
	~NuFitConfig() = default;  // destructor
	NuFitConfig(const NuFitConfig&) = default;  // copy constructor
	NuFitConfig(NuFitConfig&&) = default;  // move constructor
	NuFitConfig &operator=(const NuFitConfig&) = default;  // copy assignment
	NuFitConfig &operator=(NuFitConfig&&) = default;  // move assignment

public:  // Initialise variables to be filled by parser
	bool doToyData_ = false, doHesse = false, doMinos = false;
	double emin = 0, emax = 1e308;  // Close to numeric limit
	unsigned int nbins = 100;

	double nparams;
	std::vector<TString> param_names;
	std::vector<double> param_initial_guess, param_stepsize,
	                    param_lowerlim, param_upperlim;
};

namespace CMDLParser {
	NuFitCmdlArgs Parse(int argc, char* argv[]);
}  // namespace CMDLParser

namespace ConfigParser {
	NuFitConfig Parse(NuFitCmdlArgs);
}  // namespace ConfigParser

}  // namespace NuFitter


#endif
