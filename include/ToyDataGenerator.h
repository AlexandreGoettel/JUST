//! @file	include/ToyDataGenerator.h
//! @brief	Declaration of class to generate toy data
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef TOYDATAGENERATOR_H_
#define TOYDATAGENERATOR_H_

//============================================================================
// Includes

//============================================================================
// Class definition

namespace NuFitter {

// Forward declaration
class NuFitConfig;

class NuFitToyData {
public:
	constexpr NuFitToyData() = default;  // constructor
	~NuFitToyData() = default;  // destructor
	constexpr NuFitToyData(const NuFitToyData&) = default;  // copy constructor
	constexpr NuFitToyData(NuFitToyData&&) = default;  // move constructor
	constexpr NuFitToyData &operator=(const NuFitToyData&) = default;  // copy assignment
	constexpr NuFitToyData &operator=(NuFitToyData&&) = default;  // move assignment

public:
	bool doToyData_ = true;
};

// template <class Config>
NuFitToyData *generateToyData(const NuFitConfig *config);

}  // namespace NuFitter


#endif
