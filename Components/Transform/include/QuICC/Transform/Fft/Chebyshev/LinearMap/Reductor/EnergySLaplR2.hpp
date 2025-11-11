/**
 * @file EnergySLaplR2.hpp
 * @brief Implementation of the Chebyshev based EnergySLaplR2 reductor, with linear
 * map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGYSLAPLR2_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGYSLAPLR2_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Base/EnergySLaplR2.hpp"
#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Reductor/EnergySLaplR2viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
// #include
// "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Reductor/EnergySLaplR2viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGYSLAPLR2_HPP
