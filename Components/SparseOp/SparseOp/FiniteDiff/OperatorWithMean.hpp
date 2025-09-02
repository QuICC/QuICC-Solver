/**
 * @file OperatorWithMean.hpp
 * @brief Implementation of the generic interface for a sparse finite difference
 * operator with different treatment for the mean (l=0)
 */

#ifndef QUICC_SPARSEOP_FINITEDIFF_OPERATORWITHMEAN_HPP
#define QUICC_SPARSEOP_FINITEDIFF_OPERATORWITHMEAN_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "SparseOp/FiniteDiff/IBaseOperator.hpp"
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {
namespace SparseOp {
namespace FiniteDiff {

/// @brief Wrapper for generic Worland operator with different treatment for l=0
/// @tparam TFdBuilder builder for l!=0
/// @tparam TZeroBuilder builder for l=0
template <class TFdBuilder, class TZeroBuilder = void>
class OperatorWithMean : public IBaseOperator
{
public:
   /// @brief Pass-by-value polynomial builder ctor
   /// @param fdBuilder to be stored and used
   OperatorWithMean(TFdBuilder fdBuilder) : mFdBuilder(fdBuilder){};

   /// @brief ctor
   OperatorWithMean() = default;

   /// @brief dtor
   ~OperatorWithMean() = default;

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

template <class TFdBuilder, class TZeroBuilder>
void OperatorWithMean<TFdBuilder, TZeroBuilder>::compute(
   SparseMatrix& op, const Internal::Array& grid,
   const std::uint32_t l)
{
   assert(op.rows() == grid.size());

   if (l == 0)
   {
      if constexpr (std::is_same_v<TZeroBuilder, void>)
      {
         op.setZero();
      }
      else
      {
         TZeroBuilder fdBuilder;
         fdBuilder.compute(op, static_cast<int>(l), grid);
      }
   }
   else
   {
      mFdBuilder.compute(op, static_cast<int>(l), grid);
   }
}


} // namespace FiniteDiff
} // namespace SparseOp
} // namespace QuICC

#endif // define QUICC_SPARSEOP_FINITEDIFF_OPERATORWITHMEAN_HPP
