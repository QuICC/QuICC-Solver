/**
 * @file OperatorINWithMean.hpp
 * @brief Implementation of the generic interface for a dense Worland spectral
 * operator builder with quasi-inverse and different treatment for the mean (l=0)
 */

#ifndef QUICC_DENSESM_WORLAND_OPERATORINWITHMEAN_HPP
#define QUICC_DENSESM_WORLAND_OPERATORINWITHMEAN_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "DenseSM/IDenseSMOperator.hpp"
#include "DenseSM/Worland/GetSuperN.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {
namespace DenseSM {
/// @brief namespace for generic dense spectral Worland operator builders
namespace Worland {

/// @brief Wrapper for generic Worland operator with IN step and different
/// treatment for l=0
/// @tparam TPolyBuilder builder for l!=0
/// @tparam TINBuilder IN builder
/// @tparam TZeroBuilder builder for l=0
template <class TPolyBuilder, class TINBuilder, class TZeroBuilder = void,
   class TINZeroBuilder = void>
class OperatorINWithMean : IDenseSMOperator
{
public:
   /// @brief ctor
   OperatorINWithMean() = default;

   /// @brief dtor
   ~OperatorINWithMean() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param weights
   /// @param l
   void compute(Eigen::Ref<Matrix> op, const Internal::Array& grid,
      const Internal::Array& weights, const std::uint32_t l) final;
};

template <class TPolyBuilder, class TINBuilder, class TZeroBuilder,
   class TINZeroBuilder>
void OperatorINWithMean<TPolyBuilder, TINBuilder, TZeroBuilder,
   TINZeroBuilder>::compute(Eigen::Ref<Matrix> op, const Internal::Array& grid,
   const Internal::Array& weights, const std::uint32_t l)
{
   assert(op.rows() == grid.size());

   auto nPoly = op.cols();
   namespace ev = Polynomial::Worland::Evaluator;

   if (l == 0)
   {
      if constexpr (std::is_same_v<TZeroBuilder, void>)
      {
         static_assert(std::is_same_v<TINZeroBuilder, void>);
         op.setZero();
      }
      else
      {
         std::uint32_t extraN = getSuperN<TINZeroBuilder>();
         auto nN = nPoly + extraN;

         Internal::Matrix tOp(grid.size(), nN);
         // help compiler deducing type
         Eigen::Ref<Internal::Matrix> rTOp(tOp);

         TZeroBuilder polyBuilder;
         polyBuilder.compute(rTOp, nN, static_cast<int>(l), grid, weights,
            ev::Set());

         auto a = polyBuilder.alpha(l);
         auto b = polyBuilder.dBeta();
         TINZeroBuilder spasm(nN, nN, a, b, static_cast<int>(l),
            !static_cast<bool>(extraN));
         tOp = (spasm.mat() * tOp.transpose()).transpose();
         op = tOp.cast<MHDFloat>().leftCols(nPoly);
      }
   }
   else
   {
      std::uint32_t extraN = getSuperN<TINBuilder>();
      auto nN = nPoly + extraN;

      Internal::Matrix tOp(grid.size(), nN);
      // help compiler deducing type
      Eigen::Ref<Internal::Matrix> rTOp(tOp);

      TPolyBuilder polyBuilder;
      polyBuilder.compute(rTOp, nN, static_cast<int>(l), grid, weights,
         ev::Set());

      auto a = polyBuilder.alpha(l);
      auto b = polyBuilder.dBeta();
      TINBuilder spasm(nN, nN, a, b, static_cast<int>(l),
         !static_cast<bool>(extraN));
      tOp = (spasm.mat() * tOp.transpose()).transpose();
      op = tOp.cast<MHDFloat>().leftCols(nPoly);
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // define QUICC_DENSESM_WORLAND_OPERATORINWITHMEAN_HPP
