//! @file	src/Parser.cpp
//! @brief	Implementation of classes to parse cmdl args and config files
//! @author	Alexandre Göttel
//! @date	2021-03-18
#include <memory>
#include "Parser.h"
#include <fstream>
#include <iostream>
namespace NuFitter {

// Temporary for testing
NuFitConfig::NuFitConfig() {

}

namespace CMDLParser {

auto Parse(int argc, char* argv[]) -> NuFitCmdlArgs {
	
	auto args = std::make_unique<NuFitCmdlArgs>();

        for (int i = 1; i < argc; i++){

                if (i + 1 != argc){

                        if (strcmp(argv[i], "--config") == 0){

                                args->config1 = argv[i+1];
                                i++;
                        }
                }
        }

	return *args;
}

}  // namespace CMDLParser

namespace ConfigParser {

auto Parse(NuFitCmdlArgs args) -> NuFitConfig {
	
	std::ifstream ReadConf;
        ReadConf.open(args.config1.c_str());

        //quite rough, to be changed
        std::string a, b;
        double c, d, e, f, g, h, i, l, x;

        ReadConf >> x;
        ReadConf >> a;
        ReadConf >> b;
        ReadConf >> c;
        ReadConf >> d;
        ReadConf >> e;
        ReadConf >> f;
        ReadConf >> g;
        ReadConf >> h;
        ReadConf >> i;
        ReadConf >> l;

        std::vector<TString> s;
        s.push_back(a);
        s.push_back(b);

        std::vector<double> in;
        in.push_back(c);
        in.push_back(d);

        std::vector<double> low;
        low.push_back(e);
        low.push_back(f);

        std::vector<double> up;
        up.push_back(g);
        up.push_back(h);

        std::vector<double> step;
        step.push_back(i);
        step.push_back(l);

	auto config = std::make_unique<NuFitConfig>();
	
	config->nparams = x;
	config->param_names = s;
	config->param_initial_guess = in;
	config->param_lowerlim = low;
	config->param_upperlim = up;
	config->param_stepsize = step;

	return *config;
}

}  // namespace ConfigParser
}  // namespace NuFitter
