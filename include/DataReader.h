//! @file	include/DataReader.h
//! @brief	Declaration of classes to read data and PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef DATAREADER_H_
#define DATAREADER_H_

//============================================================================
// Includes
#include <vector>

//============================================================================
// Class definition

namespace NuFitter {

// Forward declaration
class NuFitConfig;

class NuFitData {
public:  // Constructros and operator assigments
	NuFitData() = default;  // constructor
	~NuFitData() = default;  // destructor
	NuFitData(const NuFitData&) = default;  // copy constructor
	NuFitData(NuFitData&&) = default;  // move constructor
	NuFitData &operator=(const NuFitData&) = default;  // copy assignment
	NuFitData &operator=(NuFitData&&) = default;  // move assignment

public: // Variables
	std::vector<double> data, bin_edges;
};

class NuFitPDFs {
public:  // Constructros and operator assigments
	NuFitPDFs() = default;  // constructor
	~NuFitPDFs() = default;  // destructor
	NuFitPDFs(const NuFitPDFs&) = default;  // copy constructor
	NuFitPDFs(NuFitPDFs&&) = default;  // move constructor
	NuFitPDFs &operator=(const NuFitPDFs&) = default;  // copy assignment
	NuFitPDFs &operator=(NuFitPDFs&&) = default;  // move assignment

public:  // Variables
	std::vector<double> bin_edges;
	std::vector<std::vector<double>> pdfs;
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
