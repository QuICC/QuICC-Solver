/**
 * @file Util.hpp
 * @brief Fourier operators utilities
 */
#pragma once

// External includes
//
#include <type_traits>

// Project includes
//
#include "ViewOps/Fourier/Tags.hpp"


#ifdef __CUDACC__
#define QUICC_CUDA_HOSTDEV __host__ __device__
#else
#define QUICC_CUDA_HOSTDEV
#endif

namespace QuICC {
namespace Transform {
namespace Fourier {

namespace details {

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


namespace dealias {

    /// @brief fraction of coefficients to be kept
    constexpr double rule = 2.0/3.0;

    /// @brief get starting coeff to skip copy for p bwd in place
    /// @tparam Tout
    /// @tparam Tin
    /// @tparam Order
    /// @tparam Direction
    /// @tparam Treatment
    /// @param in View
    /// @param out View
    /// @return starting mode
    template<class Tout, class Tin, std::size_t Order, class Direction, std::uint16_t Treatment>
    QUICC_CUDA_HOSTDEV inline std::size_t getMstart(Tout& out, const Tin& in)
    {
        if constexpr (std::is_same_v<Direction, bwd_t> &&
            Treatment == dealias_m &&  Order == 0)
        {
            // the diff is null and in place, we can skip the coeff to be kept
            if(out.data() == in.data())
            {
                return in.dims()[0]*rule;
            }
        }
        return 0;
    }

} // namespace dealias

} // namespace Fourier
} // namespace Transform
} // namespace QuICC
