//! @file	src/ToyDataGenerator.cpp
//! @brief	Implementation of class to generate toy data
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard includes
#include <memory>
// Project includes
#include "ToyDataGenerator.h"

namespace NuFitter{

// template <class Config>
auto generateToyData(const NuFitConfig *config) -> NuFitToyData* {
	// TODO
	auto data = std::make_unique<NuFitToyData>();
	return data.get();
}
}  // namespace NuFitter
