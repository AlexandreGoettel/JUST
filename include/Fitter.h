//! @file	include/NuFitContainer.h
//! @brief	Declaration of class to fit MC PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef NuFitContainer_H_
#define NuFitContainer_H_

//============================================================================
// Standard includes
#include <vector>
// ROOT includes
#include "TMinuit.h"
// Project includes
#include "FitResults.h"
#include "DataReader.h"
#include "ToyDataGenerator.h"

//============================================================================
// Method definitions

namespace NuFitter {

//============================================================================
// Forward declarations
class NuFitPDFs;

namespace MCFit {

class NuFitContainer {
public:  // Constructors and assigment operators
	NuFitContainer() = default;
	NuFitContainer(NuFitData*&, NuFitPDFs*&);  // constructor
	~NuFitContainer();  // destructor
	NuFitContainer(const NuFitContainer&) = delete;  // copy constructor
	NuFitContainer(NuFitContainer&&) = delete;  // move constructor
	NuFitContainer &operator=(const NuFitContainer&) = delete;  // copy assignment
	NuFitContainer &operator=(NuFitContainer&&) = delete;  // move assignment

public:
	NuFitData *data = nullptr;
	NuFitPDFs *pdfs = nullptr;
public:  // Minuit functions
	// double NLL_poisson(int, const double*);
    template <class T> T NLL_poisson(T, T);
	template <class T> T NLL_MUST(T, T);
    template <typename L> double NLL(L, int, const double*);
	std::vector<std::vector<double>> fitFunction(unsigned int, const double*);
public:  // Member variables
	unsigned int n_params, n_fixed;
	std::vector<double> efficiencies;
	std::vector<std::vector<paramData>> paramVector, paramVector_fixed;
private:  // Member variables
	std::vector<std::vector<double>> pdf_vectors, data_vector, fitValFixed;
	template <class T> bool InFitRange(T, T);
};

class MinuitManager {
public:  // Constructors and assigment operators
	MinuitManager();  // constructor
	~MinuitManager();  // destructor
	MinuitManager(const MinuitManager&) = delete;  // copy constructor
	MinuitManager(MinuitManager&&) = delete;  // move constructor
	MinuitManager &operator=(const MinuitManager&) = delete;  // copy assignment
	MinuitManager &operator=(MinuitManager&&) = delete;  // move assignment

public:  // Member variables
	TMinuit *gMinuit = nullptr;
	int errorflag, errorflag_cov;

public:  // Functions
	void initMinuit();
	void callMinuit();
	NuFitResults *getResults();
};

NuFitResults* Fit(NuFitData*&, NuFitPDFs*&);
std::vector<NuFitResults*> Fit(NuFitToyData*&, NuFitPDFs*&);
void fcn(int&, double*, double&, double*, int);

}  // namespace MCFit
}  // namespace NuFitter


#endif
