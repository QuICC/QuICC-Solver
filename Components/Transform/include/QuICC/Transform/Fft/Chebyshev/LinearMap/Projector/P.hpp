/**
 * @file P.hpp
 * @brief Implementation of the Chebyshev based P projector, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_P_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_P_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/Base/P.hpp"
#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Projector/PviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
//#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Projector/PviewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_P_HPP
