/**
 * @file Cfl.cpp
 * @brief Generic interface for global coding environment stuff
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Project includes
//
#include "Environment/Cfl.hpp"

namespace QuICC {

    double* Cfl_r=0;
	double* Cfl_dr=0;
	double* Cfl_r_ll1=0;

	int Cfl_sliceSize;
	int Cfl_nR;
	double Cfl_mcAlfvenDamping;
	double Cfl_mcAlfvenScale;
    double* newCfl_radial=0;
	double* newCfl_horizontal=0;

}
