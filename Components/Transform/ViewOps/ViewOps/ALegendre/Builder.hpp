/**
 * @file Builder.hpp
 * @brief Generic Associated Legendre operator builder
 */

#pragma once

// External includes
//
#include <Eigen/Core>

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/OuterProduct.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/ALegendre/TypeTraits.hpp"
#include "Cpp/Span.hpp"

namespace QuICC {
namespace Transform {
namespace ALegendre {

using QuICC::Patch::std::span;

/// @brief Wrapper for Eigen vector
/// @tparam T
template <class T>
using Evector = Eigen::Matrix<T, Eigen::Dynamic, 1>;


/// @brief Generic AL builder operator
/// @tparam Tview  type View of the operator
/// @tparam TPolyBuilder Polynomial builder
/// @tparam Tdata type of the computation (needed for MP)
/// @param opView operator view to be set by the builder
/// @param grid
/// @param weights
/// @param scaling
template <class Tview, class TPolyBuilder, class Tdata>
void builder(Tview opView, const Evector<Tdata>& grid,
    const Evector<Tdata>& weights, Evector<typename Tview::ScalarType>& scaling)
{
    using IndexType = typename Tview::IndexType;

    // L - harmonic degree index
    IndexType metaIdx;
    if constexpr (is_integrator_v<Tview>)
    {
        metaIdx = 1;
    }
    else if constexpr (is_projector_v<Tview>)
    {
        metaIdx = 2;
    }
    else
    {
        static_assert("builder for this type is not implemented.");
    }

    ViewBase<IndexType>& pointers = const_cast<ViewBase<IndexType>*>(opView.pointers())[metaIdx];
    ViewBase<IndexType>& indices = const_cast<ViewBase<IndexType>*>(opView.indices())[metaIdx];

    using ScalarType = typename Tview::ScalarType;
    ViewBase<ScalarType> viewData(opView.data(), opView.size());

    // Setup converters
    tempOnHostMemorySpace converterP(pointers, TransferMode::read);
    tempOnHostMemorySpace converterI(indices, TransferMode::read | TransferMode::block);
    tempOnHostMemorySpace converterD(viewData, TransferMode::write);

    // Redirect view (noop if already on cpu)
    opView = Tview(viewData.data(), viewData.size(), opView.dims(), opView.pointers(), opView.indices());

    // L - harmonic degree index
    IndexType LIdx;
    Evector<Tdata> weightsOrNot;
    if constexpr (is_integrator_v<Tview>)
    {
        LIdx = 0;
        weightsOrNot = weights;
    }
    else if constexpr (is_projector_v<Tview>)
    {
        LIdx = 1;
    }

    IndexType offSet = 0;
    IndexType layerCounter = 0;
    for(IndexType k = 0; k < opView.dims()[2]; ++k)
    {
        // check if layer is populated
        if constexpr (is_integrator_v<Tview>)
        {
            auto layerWidth = pointers[k+1] - pointers[k];
            if (layerWidth < 1)
            {
                continue;
            }
        }
        else if constexpr (is_projector_v<Tview>)
        {
            if (indices[layerCounter] != k)
            {
                continue;
            }
        }

        // temporary slice
        using slice_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
        slice_t op;

        // Build operator
        int nPoly = opView.dims()[LIdx] - k;
        op.resize(grid.size(), nPoly);
        namespace ev = Polynomial::ALegendre::Evaluator;
        // help compiler to deduce type
        Eigen::Ref<slice_t> ref(op);
        TPolyBuilder polyBuilder;
        polyBuilder.compute(ref, nPoly, static_cast<int>(k), grid, weightsOrNot, ev::Set());

        if (scaling.size() > 0)
        {
            op = op * scaling.bottomRows(op.cols()).asDiagonal();
        }

        slice_t opT;
        if constexpr (is_integrator_v<Tview>)
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

/// @brief convenience wrapper for common scalings
/// @tparam Tview  type View of the operator
/// @tparam TPolyBuilder Polynomial builder
/// @tparam Tdata type of the computation (needed for MP)
/// @tparam LlDiff order of differentiation
/// @param opView operator view to be set by the builder
/// @param grid
/// @param weights
template <class Tview, class TPolyBuilder, class Tdata, std::int16_t LlDiff = 0>
void builder(Tview opView,
    const Evector<Tdata>& grid, const Evector<Tdata>& weights)
{
    using IndexType = typename Tview::IndexType;
    // L - harmonic degree index
    IndexType LIdx;
    if constexpr (is_integrator_v<Tview>)
    {
        LIdx = 0;
    }
    else if constexpr (is_projector_v<Tview>)
    {
        LIdx = 1;
    }
    else
    {
        static_assert("builder for this type is not implemented.");
    }

    using ScalarType = typename Tview::ScalarType;
    Evector<ScalarType> scaling;
    if constexpr (LlDiff == 1)
    {
        auto L = opView.dims()[LIdx];
        scaling = Evector<ScalarType>::LinSpaced(L, 0, L-1);
        scaling = scaling.array()*(scaling.array() + 1.0);
    }
    else if constexpr (LlDiff != 0)
    {
        auto L = opView.dims()[LIdx];
        scaling = Evector<ScalarType>::LinSpaced(L, 0, L-1);
        scaling = (scaling.array()*(scaling.array() + 1.0)).pow(LlDiff);
    }

    builder<Tview, TPolyBuilder, Tdata>(opView, grid, weights, scaling);
}


} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
