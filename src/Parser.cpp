//! @file	src/Parser.cpp
//! @brief	Implementation of classes to parse cmdl args and config files
//! @author	Alexandre Göttel
//! @date	2021-03-18
#include <memory>
#include "Parser.h"

#include <iostream>
namespace NuFitter {

// Temporary for testing
NuFitConfig::NuFitConfig() {
	nparams = 2;
	param_names = {"gaus", "expo"};
	param_initial_guess = {1000, 10000};
	param_lowerlim = {0, 0};
	param_upperlim = {15000, 15000};
	param_stepsize = {15, 15};
}

namespace CMDLParser {

auto Parse(int argc, char* argv[]) -> NuFitCmdlArgs {
	// TODO
	auto args = std::make_unique<NuFitCmdlArgs>();
	return *args;
}

}  // namespace CMDLParser

namespace ConfigParser {

auto Parse(NuFitCmdlArgs args) -> NuFitConfig {
	// TODO
	auto config = std::make_unique<NuFitConfig>();
	return *config;
}

}  // namespace ConfigParser
}  // namespace NuFitter
