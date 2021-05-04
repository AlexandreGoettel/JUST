//! @file	include/Fitter.h
//! @brief	Declaration of classes to parse cmdl args and config files
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef PARSER_H_
#define PARSER_H_

//============================================================================
// Includes
#include <vector>
#include <fstream>
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

public:
	std::string gen, spec, toy;
};

class NuFitConfig {
public:
	NuFitConfig() = default;
	~NuFitConfig() = default;  // destructor
	NuFitConfig(const NuFitConfig&) = default;  // copy constructor
	NuFitConfig(NuFitConfig&&) = default;  // move constructor
	NuFitConfig &operator=(const NuFitConfig&) = default;  // copy assignment
	NuFitConfig &operator=(NuFitConfig&&) = default;  // move assignment

public:  // Initialise variables to be filled by parser
	bool doToyData_, doHesse, doMinos;
	double lifetime, mass_target, emin, emax;  // Close to numeric limit

	unsigned int nbins;
	std::string output_name;
	std::string pdf_name;
	std::string data_name;
	std::string histo_data;

	unsigned int nparams;
	unsigned int npdfs;
	std::vector<TString> param_names;
	std::vector<double> param_initial_guess, param_stepsize,
	                    param_lowerlim, param_upperlim;
};

void ErrorReading(const std::ifstream&, const std::string&);
void HelpMessage(char*);
template<typename T> void ReadAndFill_Gen(std::ifstream&, T&, T&);
int HowManySpecies(std::ifstream&, const std::string&);
template<typename T> void ReadAndFill_Spec(std::ifstream&, T&, std::vector<T>&);

namespace CMDLParser {
	NuFitCmdlArgs Parse(int argc, char* argv[]);
}  // namespace CMDLParser

namespace ConfigParser {
	NuFitConfig Parse(NuFitCmdlArgs);
}  // namespace ConfigParser

}  // namespace NuFitter


#endif
