/**
 * @file Builder.hpp
 * @brief Generic Worland operator builder
 */

#pragma once

// External includes
//
#include <Eigen/Core>

// Project includes
//
#include "DenseSM/Worland/Operator.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/Worland/TypeTraits.hpp"

namespace QuICC {
namespace Transform {
namespace Worland {


/// @brief Wrapper for Eigen vector
/// @tparam T
template <class T> using Evector = Eigen::Matrix<T, Eigen::Dynamic, 1>;


/// @brief Generic AL builder operator
/// @tparam TView  type View of the operator
/// @tparam TDenseSMBuilder Polynomial builderTDenseSMBuiler
/// @tparam TData type of the computation (needed for MP)
/// @param opView operator view to be set by the builder
/// @param grid
/// @param weights
template <class TView, class TDenseSMBuilder, class TData, class TDirection>
void builder(TView opView, const Evector<TData>& grid,
   const Evector<TData>& weights)
{
   using IndexType = typename TView::IndexType;

   // L - harmonic degree index
   IndexType metaIdx;
   if constexpr (Uniform::is_projector_v<TView, TDirection> ||
                 Uniform::is_integrator_v<TView, TDirection>)
   {
      metaIdx = 2;
   }
   else
   {
      throw std::logic_error("builder for this type is not implemented.");
   }

   ViewBase<IndexType>& pointers =
      const_cast<ViewBase<IndexType>*>(opView.pointers())[metaIdx];
   ViewBase<IndexType>& indices =
      const_cast<ViewBase<IndexType>*>(opView.indices())[metaIdx];

   using ScalarType = typename TView::ScalarType;
   ViewBase<ScalarType> viewData(opView.data(), opView.size());

   // Setup converters
   tempOnHostMemorySpace converterP(pointers, TransferMode::read);
   tempOnHostMemorySpace converterI(indices,
      TransferMode::read | TransferMode::block);
   tempOnHostMemorySpace converterD(viewData, TransferMode::write);

   // Redirect view (noop if already on cpu)
   opView = TView(viewData.data(), viewData.size(), opView.dims(),
      opView.pointers(), opView.indices());

   // L - harmonic degree index
   IndexType LIdx;
   Evector<TData> weightsOrNot;
   if constexpr (Uniform::is_integrator_v<TView, TDirection>)
   {
      LIdx = 0;
      weightsOrNot = weights;
   }
   else if constexpr (Uniform::is_projector_v<TView, TDirection>)
   {
      LIdx = 1;
   }
   else
   {
      throw std::logic_error("builder for this type is not implemented.");
   }

   IndexType offSet = 0;
   IndexType layerCounter = 0;
   for (IndexType k = 0; k < opView.dims()[2]; ++k)
   {
      // check if layer is populated
      if constexpr (Uniform::is_projector_v<TView, TDirection> ||
                    Uniform::is_integrator_v<TView, TDirection>)
      {
         if (indices[layerCounter] != k)
         {
            continue;
         }
      }

      // temporary slice
      using slice_t = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
      slice_t op;

      // Build operator
      int nPoly = opView.dims()[LIdx];
      op.resize(grid.size(), nPoly);
      // QuICC::DenseSM::Worland::Operator<ScalarType, TData, TPolyBuilder>
      // help compiler to deduce type
      Eigen::Ref<slice_t> ref(op);
      TDenseSMBuilder denseBuilder;
      denseBuilder.compute(ref, grid, weightsOrNot, k);

      slice_t opT;
      if constexpr (std::is_same_v<TDirection, fwd_t>)
      {
         opT = op.transpose();
      }
      else
      {
         opT = op;
      }

      for (int j = 0; j < opT.cols(); ++j)
      {
         for (int i = 0; i < opT.rows(); ++i)
         {
            opView(i, j, k) = opT(i, j);
         }
      }

      offSet += opT.size();

      ++layerCounter;
   }
}


} // namespace Worland
} // namespace Transform
} // namespace QuICC
