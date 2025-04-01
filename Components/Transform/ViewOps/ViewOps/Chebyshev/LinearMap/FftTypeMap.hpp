/**
 * @file FftTypeMap.hpp
 * @brief Mapping FFT backends
 */

#pragma once

// External includes
//
#include <cstdint>

// Project includes
//
#include "Fft/Fft.hpp"
#include "ViewOps/Chebyshev/LinearMap/Tags.hpp"

namespace QuICC {
namespace Transform {
namespace Chebyshev {
namespace LinearMap {

/// @brief This namespace hides implementation details
namespace details {

template <class Backend, class Tout, class Tin, class TType> struct Fft;

template <class Backend, class Tout, class Tin, class TType>
using Fft_t = typename Fft<Backend, Tout, Tin, TType>::type;

template <class Tout, class Tin, class TType>
struct Fft<viewCpu_t, Tout, Tin, TType>
{
   using type = typename QuICC::Fft::Fftw::FftOp<Tout, Tin, TType>;
};

#ifdef QUICC_HAS_CUDA_BACKEND
template <class Tout, class Tin, class TType>
struct Fft<viewGpu_t, Tout, Tin, TType>
{
   using type = typename QuICC::Fft::CuFft::FftOp<Tout, Tin, TType>;
};
#endif

#ifdef QUICC_USE_VKFFT
template <class Tout, class Tin, class TTYpe>
struct Fft<viewGpuVkFFT_t, Tout, Tin, TType>
{
   using type = typename QuICC::Fft::VkFft::FftOp<Tout, Tin, TType>;
};
#endif

} // namespace details

} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
