/**
 * @file IBaseOperator.hpp
 * @brief Implementation of the generic interface dense spectral Worland operator
 */

#ifndef QUICC_DENSEOP_WORLAND_IBASEOPERATOR_HPP
#define QUICC_DENSEOP_WORLAND_IBASEOPERATOR_HPP

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
namespace Worland {

/// @brief base class for dense Worland operator
class IBaseOperator: public IDenseOpOperator
{
public:
   /// @brief ctor
   IBaseOperator() = default;

   /// @brief dtor
   virtual ~IBaseOperator() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param weights
   /// @param l optional, operator might not depend on 3rd dimension
   virtual void compute(Eigen::Ref<Matrix> op, const Internal::Array& grid,
      const Internal::Array& weights, const std::uint32_t l = 0) = 0;
};

} // namespace Worland
} // namespace DenseOp
} // namespace QuICC

#endif // define QUICC_DENSEOP_WORLAND_IBASEOPERATOR_HPP
