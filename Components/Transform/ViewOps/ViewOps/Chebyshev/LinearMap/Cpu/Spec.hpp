/**
 * @file Spec.hpp
 * @brief Spec cpu backend
 */
#pragma once

// System includes
//
#include <cstdint>
#include <memory>

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
/// @tparam Operation type of spectral operation
/// @tparam Treatment special treatment mask, typically of mode zero
template <class Tout, class Tin, class Operation,
   std::uint16_t Treatment = 0>
class SpecOp
    : public BinaryBaseOp<SpecOp<Tout, Tin, Operation, Treatment>, Tout,
         Tin, typename Tout::ScalarType::value_type>
{
public:
   /// @brief Type of treatment mask
   static constexpr std::uint16_t TreatmentValue = Treatment;
   /// @brief Type of scale parameter, i.e. float 32/64 bits
   using ScaleType = double;
   /// @brief Constructor with user defined scaling factor
   /// @param scale
   SpecOp(ScaleType scale);
   /// @brief Default constructor
   SpecOp() = default;
   /// @brief dtor
   ~SpecOp() = default;

private:
   /// @brief Action implementation
   /// @param out output modes
   /// @param in input modes
   /// @param fftScaling fft scaling (inverse number of grid points)
   void applyImpl(Tout& out, const Tin& in, const ScaleType fftScaling);
   /// @brif Give access to base class
   friend BinaryBaseOp<SpecOp<Tout, Tin, Operation, Treatment>, Tout,
      Tin, ScaleType>;
   /// @brif Scaling factor, i.e. domain size
   ScaleType mScale{1.0};
};

} // namespace Cpu
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
