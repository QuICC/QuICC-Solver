/**
 * @file IBaseOperator.hpp
 * @brief Implementation of the generic interface dense finite difference operator
 */

#ifndef QUICC_DENSEOP_FINITEDIFF_IBASEOPERATOR_HPP
#define QUICC_DENSEOP_FINITEDIFF_IBASEOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"
#include "DenseOp/IDenseOpOperator.hpp"

namespace QuICC {
namespace DenseOp {
/// @brief namespace for generic dense Finite Differences operator builders
namespace FiniteDiff {

/// @brief base class for dense finite difference operator
class IBaseOperator: public IDenseOpOperator
{
public:
   // Operator storage type
   typedef Matrix OpType;

   /// @brief ctor
   IBaseOperator() = default;

   /// @brief dtor
   virtual ~IBaseOperator() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param l optional, operator might not depend on 3rd dimension
   virtual void compute(Eigen::Ref<Matrix> op, const Internal::Array& grid,
      const std::uint32_t l = 0) = 0;
};

} // namespace FiniteDiff
} // namespace DenseOp
} // namespace QuICC

#endif // define QUICC_DENSEOP_FINITEDIFF_IBASEOPERATOR_HPP
