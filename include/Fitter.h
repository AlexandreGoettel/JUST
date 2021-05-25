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
#include "Parser.h"

//============================================================================
// Method definitions

namespace NuFitter {

//============================================================================
// Forward declarations
class NuFitConfig;
class NuFitPDFs;

namespace MCFit {

class NuFitContainer {
public:  // Constructors and assigment operators
	NuFitContainer() = default;
	NuFitContainer(NuFitData*, NuFitPDFs*, const NuFitConfig);  // constructor
	~NuFitContainer() = default;  // destructor
	NuFitContainer(const NuFitContainer&) = default;  // copy constructor
	NuFitContainer(NuFitContainer&&) = default;  // move constructor
	NuFitContainer &operator=(const NuFitContainer&) = default;  // copy assignment
	NuFitContainer &operator=(NuFitContainer&&) = default;  // move assignment

public:
	NuFitData *data;
	NuFitPDFs *pdfs;
	NuFitConfig config;
public:  // Minuit functions
	double NLL_extended(int, const double*);
	double NLL_poisson(int, const double*);
	double fitFunction(unsigned int, unsigned int, const double*);
public:  // Member variables
	double n_params, n_fixed;
	std::vector<unsigned int> idx_map, idx_map_fixed;
	std::vector<double> efficiencies;
private:  // Member variables
	std::vector<std::vector<double>> pdf_vectors;
	std::vector<double> data_vector;
	template <class T> bool InFitRange(T, T);
};

class MinuitManager {
public:  // Constructors and assigment operators
	MinuitManager(const NuFitConfig);  // constructor
	~MinuitManager() = default;  // destructor
	MinuitManager(const MinuitManager&) = default;  // copy constructor
	MinuitManager(MinuitManager&&) = default;  // move constructor
	MinuitManager &operator=(const MinuitManager&) = default;  // copy assignment
	MinuitManager &operator=(MinuitManager&&) = default;  // move assignment

public:  // Member variables
	NuFitConfig config;
	TMinuit *gMinuit;
	int errorflag, errorflag_cov;

public:  // Functions
	void initMinuit();
	void callMinuit();
	NuFitResults getResults();
};

NuFitResults Fit(NuFitData*, NuFitPDFs*, const NuFitConfig);
NuFitResults Fit(std::vector<NuFitData*>, NuFitPDFs*,
				  const NuFitConfig);
void fcn(int&, double*, double&, double*, int);

}  // namespace MCFit
}  // namespace NuFitter


#endif
