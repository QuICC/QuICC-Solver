/**
 * @file OperatorWithMean.hpp
 * @brief Implementation of the generic interface for adense worland spectral
 * operator with different treatment for the mean (l=0)
 */

#ifndef QUICC_DENSEOP_WORLAND_OPERATORWITHMEAN_HPP
#define QUICC_DENSEOP_WORLAND_OPERATORWITHMEAN_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "DenseOp/IDenseOpOperator.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {
namespace DenseOp {
namespace Worland {

/// @brief Wrapper for generic Worland operator with different treatment for l=0
/// @tparam TPolyBuilder builder for l!=0
/// @tparam TZeroBuilder builder for l=0
template <class TPolyBuilder, class TZeroBuilder = void>
class OperatorWithMean : IDenseOpOperator
{
public:
   /// @brief Pass-by-value polynomial builder ctor
   /// @param polyBuilder to be stored and used
   OperatorWithMean(TPolyBuilder polyBuilder) : mPolyBuilder(polyBuilder){};

   /// @brief ctor
   OperatorWithMean() = default;

   /// @brief dtor
   ~OperatorWithMean() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param weights
   /// @param l
   void compute(Eigen::Ref<Matrix> op, const Internal::Array& grid,
      const Internal::Array& weights, const std::uint32_t l) final;

private:
   TPolyBuilder mPolyBuilder;
};

template <class TPolyBuilder, class TZeroBuilder>
void OperatorWithMean<TPolyBuilder, TZeroBuilder>::compute(
   Eigen::Ref<Matrix> op, const Internal::Array& grid,
   const Internal::Array& weights, const std::uint32_t l)
{
   assert(op.rows() == grid.size());

   namespace ev = Polynomial::Worland::Evaluator;
   if (l == 0)
   {
      if constexpr (std::is_same_v<TZeroBuilder, void>)
      {
         op.setZero();
      }
      else
      {
         TZeroBuilder polyBuilder;
         polyBuilder.compute(op, op.cols(), static_cast<int>(l), grid, weights,
            ev::Set());
      }
   }
   else
   {
      mPolyBuilder.compute(op, op.cols(), static_cast<int>(l), grid, weights,
         ev::Set());
   }
}


} // namespace Worland
} // namespace DenseOp
} // namespace QuICC

#endif // define QUICC_DENSEOP_WORLAND_OPERATORWITHMEAN_HPP
