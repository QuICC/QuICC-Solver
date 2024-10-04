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

/// @brief heat transport equation kernel
/// @tparam T scalar
template <class T = double> struct SphericalHeatAdvection
{
   /// @brief radial grid
   View::ViewBase<T> _grid;
   /// @brief non dimensional scaling for transport term
   T _scaling;

   /// @brief ctor
   /// @param grid
   /// @param scaling
   SphericalHeatAdvection(View::ViewBase<T> grid, T scaling) : _grid(grid), _scaling(scaling){};

   /// @brief deleted default constructor
   SphericalHeatAdvection() = delete;

   /// @brief dtor
   ~SphericalHeatAdvection() = default;

   /// @brief functor, v dot grad T minus r*v_R
   /// @param l layer index
   /// @param vR velocity, R component
   /// @param vTheta velocity, Theta component
   /// @param vPhi velocity, Phi component
   /// @param TdR temperature, R partial derivative
   /// @param TdTheta temperature, Theta partial derivative
   /// @param TdPhi temperature, Phi partial derivative
   /// @return
   QUICC_CUDA_HOSTDEV T operator()(std::uint32_t l, T vR, T vTheta, T vPhi, T TdR, T TdTheta, T TdPhi)
   {
      return _scaling * (vR*TdR + vTheta*TdTheta + vPhi*TdPhi - _grid[l]*vR);
   }
};


/// @brief Scalar R grid multiplication
/// @tparam T scalar
template <class T = double> struct MulRFunctor
{
   /// @brief radial grid
   View::ViewBase<T> _grid;
   /// @brief non dimensional scaling for transport term
   T _scaling;

   /// @brief ctor
   /// @param grid
   /// @param scaling
   MulRFunctor(View::ViewBase<T> grid, T scaling) : _grid(grid), _scaling(scaling)
   {
      assert(_grid.data() != nullptr);
   };

   /// @brief deleted default constructor
   MulRFunctor() = delete;

   /// @brief dtor
   ~MulRFunctor() = default;

   /// @brief functor, v dot grad T minus r*v_R
   /// @param l layer index
   /// @param vR velocity, R component
   /// @param vTheta velocity, Theta component
   /// @param vPhi velocity, Phi component
   /// @param TdR temperature, R partial derivative
   /// @param TdTheta temperature, Theta partial derivative
   /// @param TdPhi temperature, Phi partial derivative
   /// @return
   QUICC_CUDA_HOSTDEV T operator()(std::uint32_t l, T Rhs)
   {
      return _scaling * _grid[l] * Rhs;
   }
};



} // namespace Slicewise
} // namespace QuICC
