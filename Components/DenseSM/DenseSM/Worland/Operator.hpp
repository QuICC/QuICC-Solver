/**
 * @file Operator.hpp
 * @brief Implementation of the generic interface for a dense worland spectral
 * operator
 */

#ifndef QUICC_DENSESM_WORLAND_OPERATOR_HPP
#define QUICC_DENSESM_WORLAND_OPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "DenseSM/IDenseSMOperator.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {
namespace DenseSM {
/// @brief namespace for generic dense spectral Worland operator builders
namespace Worland {

/// @brief Wrapper for generic Worland operator
/// @tparam TPolyBuilder polynomial builder
template <class TPolyBuilder> class Operator : IDenseSMOperator
{
public:
   /// @brief ctor
   Operator() = default;

   /// @brief dtor
   ~Operator() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param weights
   /// @param l
   void compute(Eigen::Ref<Matrix> op, const Internal::Array& grid,
      const Internal::Array& weights, const std::uint32_t l) final;
};

template <class TPolyBuilder>
void Operator<TPolyBuilder>::compute(Eigen::Ref<Matrix> op,
   const Internal::Array& grid, const Internal::Array& weights,
   const std::uint32_t l)
{
   assert(op.rows() == grid.size());

   namespace ev = Polynomial::Worland::Evaluator;
   TPolyBuilder polyBuilder;
   polyBuilder.compute(op, op.cols(), static_cast<int>(l), grid, weights,
      ev::Set());
}


} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // define QUICC_DENSESM_WORLAND_OPERATOR_HPP
