/**
 * @file Functors.hpp
 * @brief Scalar functors that allow for the explicit instantiation of Cuda
 * Slicewise operators on Views.
 */
#pragma once

// System includes
//

// Project includes
//
#include "View/ViewUtils.hpp"


namespace QuICC {
/// @brief namespace for Slicewise type operations
namespace Slicewise {


/// @brief Scalar R grid multiplication
/// @tparam T scalar
template <class T = double> struct MulRFunctor
{
   /// @brief non dimensional scaling for transport term
   T _scaling;

   /// @brief ctor
   /// @param grid
   /// @param scaling
   MulRFunctor(T scaling) : _scaling(scaling) {};

   /// @brief deleted default constructor
   MulRFunctor() = delete;

   /// @brief dtor
   ~MulRFunctor() = default;

   /// @brief functor, scaled mul grid * rhs
   /// @param grid
   /// @param Rhs
   /// @return
   QUICC_CUDA_HOSTDEV T operator()(T grid, T Rhs)
   {
      return _scaling * grid * Rhs;
   }
};

/// @brief Scalar sin grid multiplication
/// @tparam T scalar
template <class T = double> struct MulSinFunctor
{
   /// @brief non dimensional scaling for transport term
   T _scaling;

   /// @brief ctor
   /// @param grid
   /// @param scaling
   MulSinFunctor(T scaling) : _scaling(scaling) {};

   /// @brief deleted default constructor
   MulSinFunctor() = delete;

   /// @brief dtor
   ~MulSinFunctor() = default;

   /// @brief functor, scaled mul grid * rhs
   /// @param grid
   /// @param Rhs
   /// @return
   QUICC_CUDA_HOSTDEV T operator()(T grid, T Rhs)
   {
      return _scaling * std::sin(grid) * Rhs;
   }
};

/// @brief Scalar cos grid multiplication
/// @tparam T scalar
template <class T = double> struct MulCosFunctor
{
   /// @brief non dimensional scaling for transport term
   T _scaling;

   /// @brief ctor
   /// @param grid
   /// @param scaling
   MulCosFunctor(T scaling) : _scaling(scaling) {};

   /// @brief deleted default constructor
   MulCosFunctor() = delete;

   /// @brief dtor
   ~MulCosFunctor() = default;

   /// @brief functor, scaled mul grid * rhs
   /// @param grid
   /// @param Rhs
   /// @return
   QUICC_CUDA_HOSTDEV T operator()(T grid, T Rhs)
   {
      return _scaling * std::cos(grid) * Rhs;
   }
};

/// @brief Scalar sin grid plus cos grid multiplication
/// @tparam T scalar
template <class T = double> struct MulSinPlusCosFunctor
{
   /// @brief non dimensional scaling for transport term
   T _scaling;

   /// @brief ctor
   /// @param grid
   /// @param scaling
   MulSinPlusCosFunctor(T scaling) : _scaling(scaling) {};

   /// @brief deleted default constructor
   MulSinPlusCosFunctor() = delete;

   /// @brief dtor
   ~MulSinPlusCosFunctor() = default;

   /// @brief functor, scaled mul grid * rhs
   /// @param grid
   /// @param Rhs
   /// @return
   QUICC_CUDA_HOSTDEV T operator()(T grid, T Lhs, T Rhs)
   {
      return _scaling * (std::sin(grid) * Lhs + std::cos(grid) * Rhs);
   }
};




} // namespace Slicewise
} // namespace QuICC
