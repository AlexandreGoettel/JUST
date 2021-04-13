//! @file	include/Fitter.h
//! @brief	Declaration of class to fit MC PDFs
//! @author	Alexandre Göttel
//! @date	2021-03-18

#ifndef FITTER_H_
#define FITTER_H_

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

class Fitter {
public:  // Constructors and assigment operators
	Fitter() = default;  // constructor
	~Fitter() = default;  // destructor
	Fitter(const Fitter&) = default;  // copy constructor
	Fitter(Fitter&&) = default;  // move constructor
	Fitter &operator=(const Fitter&) = default;  // copy assignment
	Fitter &operator=(Fitter&&) = default;  // move assignment

public:  // Functions
	// Overloaded functions to fit Data or Toy Data
	NuFitResults *Fit(NuFitData*, NuFitPDFs*, const NuFitConfig*);
	NuFitResults *Fit(std::vector<NuFitData*>, NuFitPDFs*,
		              const NuFitConfig*);

protected:  // Minuit functions
	void fcn(int&, double*, double&, double*, int);
	double NLL(int, const double*);
	double fitFunction(unsigned int, unsigned int, const double*);
	void initMinuit(const NuFitConfig*);
	void callMinuit(const NuFitConfig*);

private:  // Member variables
	std::vector<std::vector<double>> pdf_vectors;
	std::vector<double> data_vector, efficiencies;
	TMinuit *gMinuit;
	int errorflag = 0;
};

}  // namespace MCFit
}  // namespace NuFitter


#endif
