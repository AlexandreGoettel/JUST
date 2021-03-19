//! @file	include/Fitter.h
//! @brief	Declaration of classes to parse cmdl args and config files
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef PARSER_H_
#define PARSER_H_

//============================================================================
// Includes

//============================================================================
// Class definition

namespace NuFitter {
class NuFitCmdlArgs {
public:
	constexpr NuFitCmdlArgs() = default;  // constructor
	~NuFitCmdlArgs() = default;  // destructor
	constexpr NuFitCmdlArgs(const NuFitCmdlArgs&) = default;  // copy constructor
	constexpr NuFitCmdlArgs(NuFitCmdlArgs&&) = default;  // move constructor
	constexpr NuFitCmdlArgs &operator=(const NuFitCmdlArgs&) = default;  // copy assignment
	constexpr NuFitCmdlArgs &operator=(NuFitCmdlArgs&&) = default;  // move assignment
};

class NuFitConfig {
public:
	constexpr NuFitConfig() = default;  // constructor
	~NuFitConfig() = default;  // destructor
	constexpr NuFitConfig(const NuFitConfig&) = default;  // copy constructor
	constexpr NuFitConfig(NuFitConfig&&) = default;  // move constructor
	constexpr NuFitConfig &operator=(const NuFitConfig&) = default;  // copy assignment
	constexpr NuFitConfig &operator=(NuFitConfig&&) = default;  // move assignment

public:
	bool doToyData_ = true;
};

namespace CMDLParser {
	NuFitCmdlArgs Parse(int argc, char* argv[]);
}  // namespace CMDLParser

namespace ConfigParser {
	NuFitConfig *Parse(NuFitCmdlArgs);
}  // namespace ConfigParser

}  // namespace NuFitter


#endif
