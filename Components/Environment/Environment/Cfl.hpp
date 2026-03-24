/**
 * @file Cfl.hpp
 * @brief Generic interface for global coding environment stuff
 */

#ifndef QUICC_CFL_HPP
#define QUICC_CFL_HPP

// System includes
//
#include <cassert>
#include <string>
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "Environment/Typedefs.hpp"

namespace QuICC {

  extern double* Cfl_r;
	extern double* Cfl_dr;
	extern double* Cfl_r_ll1;

	extern int Cfl_sliceSize;
	extern int Cfl_nR;
	extern double Cfl_mcAlfvenDamping;
	extern double Cfl_mcAlfvenScale;
    extern double* newCfl_radial;
	extern double* newCfl_horizontal;
}

#endif // QUICC_CFL_HPP
