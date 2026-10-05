/**
 * @file Functors.hpp
 * @brief Scalar functors that allow for the explicit instantiation of Cuda
 * pointwise operators on Views.
 */
#pragma once

// System includes
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

/// @brief scalar add operation
/// specialization for complex double
#ifdef QUICC_HAS_CUDA_BACKEND
template <> struct AddFunctor<cuda::std::complex<double>>
{
   QUICC_CUDA_HOSTDEV cuda::std::complex<double> operator()(
      cuda::std::complex<double> a, cuda::std::complex<double> b)
   {
      return a + b;
   }
};
#endif

/// @brief scalar sub operation
/// @tparam T scalar
template <class T = double> struct SubFunctor
{
   QUICC_CUDA_HOSTDEV T operator()(T a, T b)
   {
      return a - b;
   }
};

/// @brief scalar add operation
/// specialization for complex double
#ifdef QUICC_HAS_CUDA_BACKEND
template <> struct SubFunctor<cuda::std::complex<double>>
{
   QUICC_CUDA_HOSTDEV cuda::std::complex<double> operator()(
      cuda::std::complex<double> a, cuda::std::complex<double> b)
   {
      return a - b;
   }
};
#endif


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

/// @tparam T scalar
template <class T = double> struct CrossCompFunctor
{
   /// @brief non dimensional scaling for transport term
   T _scaling;

   /// @brief ctor
   /// @param scaling
   CrossCompFunctor(T scaling) : _scaling(scaling){};

   /// @brief deleted default constructor
   CrossCompFunctor() = delete;

   /// @brief dtor
   ~CrossCompFunctor() = default;

   /// @brief Cross product, component wise
   /// @param uj
   /// @param uk
   /// @param vj
   /// @param vk
   /// @return i component of cross product
   QUICC_CUDA_HOSTDEV T operator()(T uj, T uk, T vj, T vk)
   {
      return _scaling * (uj * vk - uk * vj);
   }
};

/// @tparam T scalar
template <class T = double> struct ComponentDotFunctor
{
   /// @brief scaling
   T _scaling;

   /// @brief ctor
   /// @param scaling
   ComponentDotFunctor(T scaling) : _scaling(scaling){};

   /// @brief deleted default constructor
   ComponentDotFunctor() = delete;

   /// @brief dtor
   ~ComponentDotFunctor() = default;

   /// @brief Component wise dot product
   /// @param u
   /// @param v
   /// @return
   QUICC_CUDA_HOSTDEV T operator()(T u, T v)
   {
      return _scaling * (u * v);
   }

   /// @brief Component wise dot product
   /// @param u
   /// @param v
   /// @return
   QUICC_CUDA_HOSTDEV T operator()(std::complex<double> u, std::complex<double> v)
   {
      return _scaling * (u.real() * v.real() + u.imag() * v.imag());
   }
};

/// @tparam T scalar
template <class T = double> struct DotFunctor
{
   /// @brief non dimensional scaling for transport term
   T _scaling;

   /// @brief ctor
   /// @param scaling
   DotFunctor(T scaling) : _scaling(scaling){};

   /// @brief deleted default constructor
   DotFunctor() = delete;

   /// @brief dtor
   ~DotFunctor() = default;

   /// @brief Dot product
   /// @param ui
   /// @param uj
   /// @param uk
   /// @param vi
   /// @param vj
   /// @param vk
   /// @return
   QUICC_CUDA_HOSTDEV T operator()(T ui, T uj, T uk, T vi, T vj, T vk)
   {
      return _scaling * (ui * vi + uj * vj + uk * vk);
   }
};


} // namespace Pointwise
} // namespace QuICC
