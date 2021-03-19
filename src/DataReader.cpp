//! @file	src/DataReader.cpp
//! @brief	Implementation of classes to read data and PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard includes
#include <memory>
// Project includes
#include "DataReader.h"


namespace NuFitter {
namespace Data {

auto Read(const NuFitConfig *config) -> NuFitData* {
	// TODO
	auto data = std::make_unique<NuFitData>();
	return data.get();
}

}  // namespace Data

namespace PDFs {

auto Read(const NuFitConfig *config) -> NuFitPDFs* {
	// TODO
	auto data = std::make_unique<NuFitPDFs>();
	return data.get();
}

}  // namespace PDFs
}  // namespace NuFitter
