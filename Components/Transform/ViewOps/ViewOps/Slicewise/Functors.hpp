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



} // namespace Slicewise
} // namespace QuICC
