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
/// @tparam Treatment special treatment mask
template <class Tout, class Tin, class Operation, std::uint16_t Treatment>
class SpecOp : public BinaryBaseOp<SpecOp<Tout, Tin, Operation, Treatment>,
                  Tout, Tin, typename Tout::ScalarType::value_type>
{
public:
   /// @brief Type of treatment mask
   static constexpr std::uint16_t TreatmentValue = Treatment;
   /// @brief Type of scale parameter, i.e. float 32/64 bits
   using ScaleType = double;
   /// @brief Constructor with bounds and user defined scaling factor
   /// @param lower  lower bound
   /// @param upper  upper bound
   /// @param scale
   SpecOp(const double lower, const double upper, ScaleType scale);
   /// @brief constructor with bounds
   /// @param lower  lower bound
   /// @param upper  upper bound
   SpecOp(const double lower, const double upper);
   /// @brief dtor
   ~SpecOp() = default;

private:
   /// @brief Multiply by Y
   std::size_t multiplyByY(typename Tout::ScalarType* const out,
      typename Tout::ScalarType* in, const std::size_t Nout,
      const std::size_t Nin, const double c, const std::size_t shiftIn = 0);

   /// @brief Differentiate
   void differentiate(typename Tout::ScalarType* const out,
      const std::size_t Nout);

   /// @brief Action implementation
   /// @param out output modes
   /// @param in input modes
   /// @param fftScaling fft scaling (inverse number of grid points)
   void applyImpl(Tout& out, const Tin& in, const ScaleType fftScaling);
   /// @brif Give access to base class
   friend BinaryBaseOp<SpecOp<Tout, Tin, Operation, Treatment>, Tout, Tin,
      ScaleType>;
   /// @brif Lower bound
   double mLower;
   /// @brif Upper bound
   double mUpper;
   /// @brif Scaling factor, i.e. domain size
   ScaleType mScale{1.0};
};

} // namespace Cpu
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
