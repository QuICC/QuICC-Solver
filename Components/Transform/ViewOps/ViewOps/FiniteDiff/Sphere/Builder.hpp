/**
 * @file Builder.hpp
 * @brief Generic Finite Differences operator builder
 */

#pragma once

// System includes
//
#include <Eigen/Core>

// Project includes
//
#include "SparseOp/FiniteDiff/Operator.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/FiniteDiff/Sphere/TypeTraits.hpp"

namespace QuICC {
namespace Transform {
namespace FiniteDiff {
namespace Sphere {

/// @brief Generic Finite Difference builder operator
/// @tparam TView  type View of the operator
/// @tparam TSparseOpBuilder FD builder
/// @tparam TDirection
template <class TView, class TSparseOpBuilder, class TDirection> class Builder
{
public:
   /// @brief Pass-by-value dense builder ctor
   /// @param denseBuilder to be stored and used
   Builder(TSparseOpBuilder denseBuilder) : mSparseBuilder(denseBuilder){};

   /// @brief default ctor
   Builder() = default;

   /// @brief dtor
   ~Builder() = default;

   void compute(TView opView, const Internal::Array& grid);

private:
   TSparseOpBuilder mSparseBuilder;
};


template <class TView, class TSparseOpBuilder, class TDirection>
void Builder<TView, TSparseOpBuilder, TDirection>::compute(TView opView,
   const Internal::Array& grid)
{
   using IndexType = typename TView::IndexType;

   // L - harmonic degree index
   IndexType metaIdx;
   /// \todo direction might not be needed
   if constexpr (is_projector_v<TView, TDirection> ||
                 is_integrator_v<TView, TDirection>)
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
   if constexpr (is_integrator_v<TView, TDirection>)
   {
      LIdx = 0;
   }
   else if constexpr (is_projector_v<TView, TDirection>)
   {
      LIdx = 1;
   }
   else
   {
      throw std::logic_error("builder for this type is not implemented.");
   }

   IndexType layerCounter = 0;
   for (IndexType k = 0; k < opView.dims()[2]; ++k)
   {
      // check if layer is populated
      if constexpr (is_projector_v<TView, TDirection> ||
                    is_integrator_v<TView, TDirection>)
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
      using sparse_t = Eigen::SparseMatrix<ScalarType>;
      sparse_t op;

      // Build operator
      op.resize(grid.size(), grid.size());
      // QuICC::SparseOp::FiniteDiff::Operator<ScalarType, TData, TFdBuilder>
      // help compiler to deduce type
      mSparseBuilder.compute(op, grid, k);

      slice_t opT;
      opT = op;

      for (int j = 0; j < opT.cols(); ++j)
      {
         for (int i = 0; i < opT.rows(); ++i)
         {
            opView(i, j, k) = opT(i, j);
         }
      }

      ++layerCounter;
   }
}


} // namespace Sphere
} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC
