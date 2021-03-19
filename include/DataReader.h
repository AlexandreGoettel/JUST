//! @file	include/DataReader.h
//! @brief	Declaration of classes to read data and PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef DATAREADER_H_
#define DATAREADER_H_

//============================================================================
// Includes
// #include "Parser.h"

//============================================================================
// Class definition

namespace NuFitter {

// Forward declaration
class NuFitConfig;

class NuFitData {
public:
	constexpr NuFitData() = default;  // constructor
	~NuFitData() = default;  // destructor
	constexpr NuFitData(const NuFitData&) = default;  // copy constructor
	constexpr NuFitData(NuFitData&&) = default;  // move constructor
	constexpr NuFitData &operator=(const NuFitData&) = default;  // copy assignment
	constexpr NuFitData &operator=(NuFitData&&) = default;  // move assignment
};

class NuFitPDFs {
public:
	constexpr NuFitPDFs() = default;  // constructor
	~NuFitPDFs() = default;  // destructor
	constexpr NuFitPDFs(const NuFitPDFs&) = default;  // copy constructor
	constexpr NuFitPDFs(NuFitPDFs&&) = default;  // move constructor
	constexpr NuFitPDFs &operator=(const NuFitPDFs&) = default;  // copy assignment
	constexpr NuFitPDFs &operator=(NuFitPDFs&&) = default;  // move assignment
};

namespace Data {
	// template <class Config>
	NuFitData *Read(const NuFitConfig *config);
}  // namespace CMDLParser

namespace PDFs {
	// template <class Config>
	NuFitPDFs *Read(const NuFitConfig *config);
}  // namespace ConfigParser

}  // namespace NuFitter


#endif
