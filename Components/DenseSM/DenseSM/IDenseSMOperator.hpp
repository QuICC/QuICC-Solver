/**
 * @file IDenseSMOperator.hpp
 * @brief Implementation of the generic interface dense spectral operator
 */

#ifndef QUICC_DENSESM_IDENSESMOPERATOR_HPP
#define QUICC_DENSESM_IDENSESMOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//

namespace QuICC {
/// @brief namespace for generic dense spectral operator builders
namespace DenseSM {

/// @brief base class for dense operator
class IDenseSMOperator
{
public:
   /// @brief ctor
   IDenseSMOperator() = default;

   /// @brief dtor
   virtual ~IDenseSMOperator() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param weights
   /// @param l optional, operator might not depend on 3rd dimension
   virtual void compute(Eigen::Ref<Matrix> op, const Internal::Array& grid,
      const Internal::Array& weights, const std::uint32_t l = 0) = 0;
};

} // namespace DenseSM
} // namespace QuICC

#endif // define QUICC_DENSESM_IDENSESMOPERATOR_HPP
