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
#include "ViewOps/Fourier/Tags.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {

/// @brief This namespace hides implementation details
namespace details {

template <class Backend, class Tout, class Tin> struct Fft;

template <class Backend, class Tout, class Tin>
using Fft_t = typename Fft<Backend, Tout, Tin>::type;

template <class Tout, class Tin> struct Fft<viewCpu_t, Tout, Tin>
{
   using type = typename QuICC::Fft::Fftw::FftOp<Tout, Tin>;
};

#ifdef QUICC_HAS_CUDA_BACKEND
template <class Tout, class Tin> struct Fft<viewGpu_t, Tout, Tin>
{
   using type = typename QuICC::Fft::CuFft::FftOp<Tout, Tin>;
};
#endif

#ifdef QUICC_USE_VKFFT
template <class Tout, class Tin> struct Fft<viewGpuVkFFT_t, Tout, Tin>
{
   using type = typename QuICC::Fft::VkFft::FftOp<Tout, Tin>;
};
#endif

} // namespace details

} // namespace Fourier
} // namespace Transform
} // namespace QuICC
