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
#include <cstring>
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
	std::string gen, spec, toy, output_name;
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
	bool doHesse = false, doMinos = false;
	unsigned int ToyData = 0; //0 by default means no toyfit
	double emin = 0., emax = 5000.;  // Close to numeric limit
	std::string likelihood = "poisson";

	// Parameters with no standard value must be given in the config file
	double lifetime, densityLS, radius, mass_target, daq_time, exposure;  // mass_target in kt
	std::string pdf_name, data_name, output_name;

	std::vector<std::string> data_hist_names;
	std::vector<unsigned int> nbins;

	unsigned int npdfs, nparams;
	std::vector<TString> param_names, pdf_names;
	std::vector<double> param_initial_guess, param_stepsize,
	                    param_lowerlim, param_upperlim, param_eff;
	std::vector<unsigned int> param_fixed, hist_id, nSp_histos;

	// Toy_rates config file
	unsigned int npdfs_toy;
	std::vector<TString> pdf_names_toy;
	std::vector<double> param_initial_guess_toy, param_eff_toy;
	std::vector<unsigned int> hist_id_toy, nSp_histos_toy;

};

void ErrorReading(std::string&);
void HelpMessage(char*);
template<typename T> void ReadAndFill_Spec(std::ifstream&, T&, std::vector<T>&);
inline std::string GetValue(std::string&, std::string&);

namespace CMDLParser {
	NuFitCmdlArgs Parse(int argc, char* argv[]);
}  // namespace CMDLParser

namespace ConfigParser {
	void ParseGenOpts(std::unique_ptr<NuFitConfig>&, std::string);
	void ParseSpeciesList(std::unique_ptr<NuFitConfig>&, std::string);
	void ParseToyRates(std::unique_ptr<NuFitConfig>&, std::string);
	NuFitConfig Parse(NuFitCmdlArgs);
}  // namespace ConfigParser

}  // namespace NuFitter


#endif
