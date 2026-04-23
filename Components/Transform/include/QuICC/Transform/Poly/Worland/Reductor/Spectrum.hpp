/**
 * @file Spectrum.hpp
 * @brief Implementation of the Worland based R^2 power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_SPECTRUM_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_SPECTRUM_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/Spectrum.hpp"
#include "QuICC/Transform/Wrappers/Worland/Reductor/SpectrumviewCpu_t.hpp.inc"
//#ifdef QUICC_USE_KOKKOS
//#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/Spectrum.hpp"
//#endif
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Reductor/SpectrumviewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_SPECTRUM_HPP
