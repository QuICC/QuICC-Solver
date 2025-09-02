/**
 * @file IBaseOperator.hpp
 * @brief Implementation of the generic interface sparse finite difference operator
 */

#ifndef QUICC_SPARSEOP_FINITEDIFF_IBASEOPERATOR_HPP
#define QUICC_SPARSEOP_FINITEDIFF_IBASEOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"
#include "SparseOp/ISparseOpOperator.hpp"

namespace QuICC {
/// @brief namespace for generic sparse spectral operator builders
namespace SparseOp {
namespace FiniteDiff {

/// @brief base class for sparse finite difference operator
class IBaseOperator: public ISparseOpOperator
{
public:
   // Operator storage type
   typedef SparseMatrix OpType;

   /// @brief ctor
   IBaseOperator() = default;

   /// @brief dtor
   virtual ~IBaseOperator() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param l optional, operator might not depend on 3rd dimension
   virtual void compute(SparseMatrix& op, const Internal::Array& grid,
      const std::uint32_t l = 0) = 0;
};

} // namespace FiniteDiff
} // namespace SparseOp
} // namespace QuICC

#endif // define QUICC_SPARSEOP_FINITEDIFF_IBASEOPERATOR_HPP
