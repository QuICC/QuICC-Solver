/**
 * @file Operator.hpp
 * @brief Implementation of the generic interface for a sparse finite difference operator
 */

#ifndef QUICC_SPARSEOP_FINITEDIFF_OPERATOR_HPP
#define QUICC_SPARSEOP_FINITEDIFF_OPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "SparseOp/ISparseOpOperator.hpp"
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {
namespace SparseOp {
/// @brief namespace for generic sparse Finite Differences operator builders
namespace FiniteDiff {

/// @brief Wrapper for generic Finite Differences operator
/// @tparam TFdBuilder finite difference builder
template <class TFdBuilder> class Operator : ISparseOpOperator
{
public:
   /// @brief Pass-by-value polynomial builder ctor
   /// @param fdBuilder to be stored and used
   Operator(TFdBuilder fdBuilder) : mFdBuilder(fdBuilder){};

   /// @brief default ctor
   Operator() = default;

   /// @brief dtor
   ~Operator() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param weights
   /// @param l
   void compute(SparseMatrix& op, const Internal::Array& grid,
      const std::uint32_t l) final;

private:
   TFdBuilder mFdBuilder;
};

template <class TFdBuilder>
void Operator<TFdBuilder>::compute(SparseMatrix& op,
   const Internal::Array& grid,
   const std::uint32_t l)
{
   assert(op.rows() == grid.size());

   mFdBuilder.compute(op, static_cast<int>(l), grid);
}


} // namespace FiniteDiff
} // namespace SparseOp
} // namespace QuICC

#endif // define QUICC_SPARSEOP_FINITEDIFF_OPERATOR_HPP
