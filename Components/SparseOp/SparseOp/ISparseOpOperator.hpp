/**
 * @file ISparseOpOperator.hpp
 * @brief Implementation of the generic interface sparse operator
 */

#ifndef QUICC_SPARSEOP_ISPARSEOPOPERATOR_HPP
#define QUICC_SPARSEOP_ISPARSEOPOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {
/// @brief namespace for generic sparse spectral operator builders
namespace SparseOp {

/// @brief base class for sparse operator
class ISparseOpOperator
{
public:
   /// @brief ctor
   ISparseOpOperator() = default;

   /// @brief dtor
   virtual ~ISparseOpOperator() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param l optional, operator might not depend on 3rd dimension
   virtual void compute(SparseMatrix& op, const Internal::Array& grid,
      const std::uint32_t l = 0) = 0;
};

} // namespace SparseOp
} // namespace QuICC

#endif // define QUICC_SPARSEOP_ISPARSEOPOPERATOR_HPP
