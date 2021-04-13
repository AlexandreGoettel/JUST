//! @file	include/NuFitResults.h
//! @brief	Declaration of class to contain fit results
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef NuFitResults_H_
#define NuFitResults_H_

//============================================================================
// Includes

//============================================================================
// Class definition

namespace NuFitter {

class NuFitResults {
// Constructors and operator assignments
public:
	constexpr NuFitResults() = default;  // constructor
	~NuFitResults() = default;  // destructor
	constexpr NuFitResults(const NuFitResults&) = default;  // copy constructor
	constexpr NuFitResults(NuFitResults&&) = default;  // move constructor
	constexpr NuFitResults &operator=(const NuFitResults&) = default;  // copy assignment
	constexpr NuFitResults &operator=(NuFitResults&&) = default;  // move assignment

// Functions
public:
	void addResults(NuFitResults*);

};

}  // namespace NuFitter


#endif
