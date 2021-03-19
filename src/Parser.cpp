//! @file	src/Parser.cpp
//! @brief	Implementation of classes to parse cmdl args and config files
//! @author	Alexandre Göttel
//! @date	2021-03-18
#include <memory>
#include "Parser.h"

namespace NuFitter {
namespace CMDLParser {

auto Parse(int argc, char* argv[]) -> NuFitCmdlArgs {
	// TODO
	auto args = std::make_unique<NuFitCmdlArgs>();
	return *args;
}

}  // namespace CMDLParser

namespace ConfigParser {

auto Parse(NuFitCmdlArgs args) -> NuFitConfig* {
	// TODO
	auto config = std::make_unique<NuFitConfig>();
	return config.get();
}

}  // namespace ConfigParser
}  // namespace NuFitter
