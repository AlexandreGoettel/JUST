//! @file	include/Fitter.cpp
//! @brief	Implementation of class to fit MC PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

// Standard includes
#include <memory>
// Project includes
#include "Fitter.h"


namespace NuFitter {
namespace MCFit {

auto Fit(NuFitData *data, NuFitPDFs *pdfs, const NuFitConfig *config) -> NuFitResults* {
	// TODO
	auto results = std::make_unique<NuFitResults>();
	return results.get();
}

auto Fit(NuFitToyData *data, NuFitPDFs *pdfs, const NuFitConfig *config) -> NuFitResults* {
	// TODO
	auto results = std::make_unique<NuFitResults>();
	return results.get();
}

}  // namespace MCFit
}  // namespace NuFitter
