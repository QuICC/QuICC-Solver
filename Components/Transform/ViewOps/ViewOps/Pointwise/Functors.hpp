/**
 * @file Functors.hpp
 * @brief Scalar functors that allow for the explicit instantiation of Cuda
 * pointwise operators on Views.
 */
#pragma once

// External includes
//
#include <complex>
#ifdef QUICC_HAS_CUDA_BACKEND
#include <cuda/std/complex>
#endif

// Project includes
//
#include "View/ViewUtils.hpp"


namespace QuICC {
/// @brief namespace for Pointwise type operations
namespace Pointwise {

/// @brief scalar add operation
/// @tparam T scalar
template <class T = double> struct AddFunctor
{
   QUICC_CUDA_HOSTDEV T operator()(T a, T b)
   {
      return a + b;
   }
};


/// @brief scalar square operation
/// @tparam T scalar
template <class T = double> struct SquareFunctor
{
   QUICC_CUDA_HOSTDEV T operator()(T in)
   {
      return in * in;
   }
};

/// @brief scalar square absolute value operation
/// @tparam T scalar
template <class T = double> struct Abs2Functor
{
   QUICC_CUDA_HOSTDEV T operator()(std::complex<T> in)
   {
#ifdef QUICC_HAS_CUDA_BACKEND
      cuda::std::complex<T>* ptr =
         reinterpret_cast<cuda::std::complex<T>*>(&in);
      return (*ptr).real() * (*ptr).real() + (*ptr).imag() * (*ptr).imag();
#else
      return in.real() * in.real() + in.imag() * in.imag();
#endif
   }
};

} // namespace Pointwise
} // namespace QuICC
