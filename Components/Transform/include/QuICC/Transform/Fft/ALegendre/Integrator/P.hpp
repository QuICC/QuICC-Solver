/**
 * @file P.hpp
 * @brief Implementation of the Worland based P integrator
 */

#ifndef QUICC_TRANSFORM_FFT_ALEGENDRE_INTEGRATOR_P_HPP
#define QUICC_TRANSFORM_FFT_ALEGENDRE_INTEGRATOR_P_HPP

// External includes
//

// Project includes
//

#include "QuICC/Transform/Fft/ALegendre/Integrator/Base/P.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/PviewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_P_HPP