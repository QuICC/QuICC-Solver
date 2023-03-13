/**
 * @file Util.hpp
 * @brief Fourier operators utilities
 */
#pragma once

// External includes
//

// Project includes
//

#include <type_traits>

#ifdef __CUDACC__
#define QUICC_CUDA_HOSTDEV __host__ __device__
#else
#define QUICC_CUDA_HOSTDEV
#endif

namespace QuICC {
namespace Transform {
namespace Fourier {

namespace details
{
    /// Generic
    template <std::size_t n>
    QUICC_CUDA_HOSTDEV inline double fast_pow(const double i);

    template <>
    QUICC_CUDA_HOSTDEV inline double fast_pow<0>(const double i)
    {
        return 1;
    }

    template <>
    QUICC_CUDA_HOSTDEV inline double fast_pow<1>(const double i)
    {
        return i;
    }

    template <>
    QUICC_CUDA_HOSTDEV inline double fast_pow<2>(const double i)
    {
        return i*i;
    }

    template <>
    QUICC_CUDA_HOSTDEV inline double fast_pow<3>(const double i)
    {
        return i*i*i;
    }

    template <>
    QUICC_CUDA_HOSTDEV inline double fast_pow<4>(const double i)
    {
        auto i2 = i*i;
        return i2*i2;
    }
} // namespace details

} // namespace Fourier
} // namespace Transform
} // namespace QuICC
