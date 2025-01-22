/**
 * @file IDenseOpOperator.hpp
 * @brief Implementation of the generic interface dense spectral operator
 */

#ifndef QUICC_DENSEOP_IDENSEOPOPERATOR_HPP
#define QUICC_DENSEOP_IDENSEOPOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//

namespace QuICC {
/// @brief namespace for generic dense spectral operator builders
namespace DenseOp {

/// @brief base class for dense operator
class IDenseOpOperator
{
public:
   /// @brief ctor
   IDenseOpOperator() = default;

   /// @brief dtor
   virtual ~IDenseOpOperator() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param weights
   /// @param l optional, operator might not depend on 3rd dimension
   virtual void compute(Eigen::Ref<Matrix> op, const Internal::Array& grid,
      const Internal::Array& weights, const std::uint32_t l = 0) = 0;
};

} // namespace DenseOp
} // namespace QuICC

#endif // define QUICC_DENSEOP_IDENSEOPOPERATOR_HPP
