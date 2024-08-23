/**
 * @file Builder.hpp
 * @brief Generic Worland operator builder
 */

#pragma once

// System includes
//
#include <Eigen/Core>

// Project includes
//
#include "DenseSM/Worland/Operator.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/Worland/TypeTraits.hpp"

namespace QuICC {
namespace Transform {
namespace Worland {

/// @brief Generic JW builder operator
/// @tparam TView  type View of the operator
/// @tparam TDenseSMBuilder Polynomial builder
/// @tparam TDirection
template <class TView, class TDenseSMBuilder, class TDirection> class Builder
{
public:
   /// @brief Pass-by-value dense builder ctor
   /// @param denseBuilder to be stored and used
   Builder(TDenseSMBuilder denseBuilder) : mDenseBuilder(denseBuilder){};

   /// @brief default ctor
   Builder() = default;

   /// @brief dtor
   ~Builder() = default;

   void compute(TView opView, const Internal::Array& grid,
      const Internal::Array& weights);

private:
   TDenseSMBuilder mDenseBuilder;
};


template <class TView, class TDenseSMBuilder, class TDirection>
void Builder<TView, TDenseSMBuilder, TDirection>::compute(TView opView,
   const Internal::Array& grid, const Internal::Array& weights)
{
   using IndexType = typename TView::IndexType;

   // L - harmonic degree index
   IndexType metaIdx;
   /// \todo direction might not be needed, check after triangular truncation is
   /// implemented
   if constexpr (Uniform::is_projector_v<TView, TDirection> ||
                 Uniform::is_integrator_v<TView, TDirection>)
   {
      metaIdx = 2;
   }
   else
   {
      throw std::logic_error("builder for this type is not implemented.");
   }

   using namespace QuICC::Memory;
   using namespace QuICC::View;

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
   if constexpr (Uniform::is_integrator_v<TView, TDirection>)
   {
      LIdx = 0;
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
         if (layerCounter >= indices.size())
         {
            break;
         }
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
      mDenseBuilder.compute(ref, grid, weights, k);

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
