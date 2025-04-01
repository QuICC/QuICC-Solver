/**
 * @file Grid.hpp
 * @brief Grid cpu backend
 */
#pragma once

// System includes
//
#include <cstdint>
#include <memory>
#include <vector>

// Project includes
//
#include "Operator/Binary.hpp"

namespace QuICC {
namespace Transform {
/// @brief namespace for Chebyshev based operators
namespace Chebyshev {
/// @brief namespace for LinearMap Chebyshev operators (y = ax + b)
namespace LinearMap {
/// @brief Cpu backend namespace
namespace Cpu {

using namespace QuICC::Operator;

/// @brief Derived classes implement the spectral operations in modal space
/// the padded region is set to zero
/// @tparam Tout output modes type
/// @tparam Tin input modes type
/// @tparam Operation Type of grid operation
/// @tparam Treatment special treatment mask, typically of mode zero
template <class Tout, class Tin, class Operation,
   std::uint16_t Treatment = 0>
class GridOp
    : public BinaryBaseOp<GridOp<Tout, Tin, Operation, Treatment>, Tout,
         Tin, typename Tout::ScalarType::value_type>
{
public:
   /// @brief Type of treatment mask
   static constexpr std::uint16_t TreatmentValue = Treatment;
   /// @brief Type of scale parameter, i.e. float 32/64 bits
   using ScaleType = double;
   /// @brief Constructor with bounds user defined scaling factor
   /// @param lower  lower bound
   /// @param upper  upper bound
   /// @param scale
   GridOp(const double lower, const double upper, ScaleType scale);
   /// @brief constructor with bounds
   /// @param lower  lower bound
   /// @param upper  upper bound
   GridOp(const double lower, const double upper);
   /// @brief dtor
   ~GridOp() = default;

private:
   /// @brief Action implementation
   /// @param out output modes
   /// @param in input modes
   /// @param fftScaling fft scaling (inverse number of grid points)
   void applyImpl(Tout& out, const Tin& in, const ScaleType fftScaling);
   /// @brif Give access to base class
   friend BinaryBaseOp<GridOp<Tout, Tin, Operation, Treatment>, Tout,
      Tin, ScaleType>;
   /// @brief Lower bound
   double mLower;
   /// @brief Upper bound
   double mUpper;
   /// @brief
   std::vector<double> mGridScaler;
   /// @brief Scaling factor, i.e. domain size
   ScaleType mScale{1.0};
};

} // namespace Cpu
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
