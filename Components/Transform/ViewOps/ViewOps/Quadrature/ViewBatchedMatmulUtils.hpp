/**
 * @file ViewBatchedMatmulUtils.hpp
 * @brief Utilities to get sizes of batched matmul for different View layouts
 */
#pragma once

// External includes
//

// Project includes
//
#include "Std/Cuda/Utility.hpp"
#include "Std/Span.hpp"
#include "View/View.hpp"

namespace QuICC {
namespace Transform {
namespace Quadrature {

/// @brief POD to return op metadata as a single object
/// @tparam T
template <class T> struct opMeta
{
   std::size_t dataSize;
   T pointersSize, indicesSize, idx;
};

/// @brief compute operator metadata based on type, dimensions and number of
/// layers
/// @tparam Top
/// @param dimensions
/// @param nLayers
/// @return
template <class Top>
QUICC_CUDA_HOSTDEV inline opMeta<typename Top::IndexType> getOpMeta(
   span<const typename Top::IndexType> dimensions,
   span<const typename Top::IndexType> layers)
{
   using namespace QuICC::View;
   auto nLayers = layers.size();
   using IndexType = typename Top::IndexType;
   opMeta<IndexType> meta;
   using opLevelType = typename Top::LevelType;
   if constexpr (std::is_same_v<opLevelType, CSL3D::level>)
   {
      // uniform truncation projector/integrator
      // dim 0 - Nr - radial points/modes
      // dim 1 - R  - radial modes/points
      // dim 2 - L  - harmonic degree
      auto M = dimensions[0];
      auto K = dimensions[1];
      meta.pointersSize = 2;
      meta.indicesSize = nLayers;
      meta.idx = 2;
      meta.dataSize = M * K * nLayers;
   }
   // else if constexpr (is_integrator_v<Top>)
   else if constexpr (std::is_same_v<opLevelType, S1CLCSC3D::level> ||
                      std::is_same_v<opLevelType, S1CLCSC3DJIK::level>)
   {
      /// dim 0 - L  - harmonic degree
      /// dim 1 - Nl - longitudinal points
      /// dim 2 - M  - harmonic order
      auto K = dimensions[0];
      auto M = dimensions[1];
      auto P = dimensions[2];
      meta.pointersSize = P + 1;
      meta.indicesSize = nLayers * M;
      meta.idx = 1;

      std::size_t varSize = 0;
      for (IndexType i = 0; i < nLayers; ++i)
      {
         auto nPoly = K - layers[i];
         varSize += nPoly;
      }
      meta.dataSize = M * varSize;
   }
   // else if constexpr (ALegendre::is_integrator_v<Top>)
   else if constexpr (std::is_same_v<opLevelType, CS1RL3D::level> ||
                      std::is_same_v<opLevelType, CS1RL3DJIK::level>)
   {
      /// dim 0 - Nl - longitudinal points
      /// dim 1 - L  - harmonic degree
      /// dim 2 - M  - harmonic order
      auto M = dimensions[0];
      auto K = dimensions[1];
      meta.pointersSize = 2;
      meta.indicesSize = nLayers;
      meta.idx = 2;

      std::size_t varSize = 0;
      for (IndexType i = 0; i < nLayers; ++i)
      {
         auto nPoly = K - layers[i];
         varSize += nPoly;
      }
      meta.dataSize = M * varSize;
   }
   else
   {
      static_assert(std::is_same_v<opLevelType, void>,
         "ctor for this type is not implemented yet");
   }

   return meta;
}

template <class Top>
QUICC_CUDA_HOSTDEV inline void setIndicesAndPointers(
   View::ViewBase<typename Top::IndexType> pointers[],
   View::ViewBase<typename Top::IndexType> indices[],
   span<const typename Top::IndexType> dimensions,
   span<const typename Top::IndexType> layers)
{
   using namespace QuICC::View;
   using IndexType = typename Top::IndexType;
   using opLevelType = typename Top::LevelType;
   auto nLayers = layers.size();

   // Set up pointers / indices for operator
   if constexpr (std::is_same_v<opLevelType, CSL3D::level>)
   {
      IndexType metaIdx = 2;
      // uniform truncation projector/integrator
      pointers[metaIdx][0] = 0;
      pointers[metaIdx][1] = indices[metaIdx].size();

      // Full slice per layer
      for (IndexType p = 0; p < nLayers; ++p)
      {
         // set index
         indices[metaIdx][p] = layers[p];
      }
   }
   // else if constexpr (ALegendre::is_integrator_v<Top>)
   else if constexpr (std::is_same_v<opLevelType, S1CLCSC3D::level> ||
                      std::is_same_v<opLevelType, S1CLCSC3DJIK::level>)
   {
      /// dim 0 - L  - harmonic degree
      /// dim 1 - Nl - longitudinal points
      /// dim 2 - M  - harmonic order
      auto M = dimensions[1];
      auto P = dimensions[2];

      IndexType metaIdx = 1;
      pointers[metaIdx][0] = 0;
      IndexType nCols = 0;
      IndexType layerCtr = 0;
      for (IndexType i = 0; i < P; ++i)
      {
         IndexType blockWidth = 0;
         if (static_cast<IndexType>(layers[layerCtr] == i))
         {
            blockWidth = M;
            ++layerCtr;
         }

         for (IndexType idx = 0; idx < blockWidth; idx++)
         {
            indices[metaIdx][nCols + idx] = idx;
         }

         nCols += blockWidth;
         pointers[metaIdx][i + 1] = nCols;
      }
   }
   // else if constexpr (ALegendre::is_projector_v<Top>)
   else if constexpr (std::is_same_v<opLevelType, CS1RL3D::level> ||
                      std::is_same_v<opLevelType, CS1RL3DJIK::level>)
   {
      IndexType metaIdx = 2;
      pointers[metaIdx][0] = 0;
      pointers[metaIdx][1] = indices[metaIdx].size();

      // Copy operators matrix by matrix to view
      for (IndexType p = 0; p < nLayers; ++p)
      {
         // set index
         indices[metaIdx][p] = layers[p];
      }
   }
   else
   {
      static_assert(std::is_same_v<opLevelType, void>,
         "ctor for this type is not implemented yet");
   }
}

/// @brief Get the pointers metadata for the modes View.
/// From this we can check if the layer is populated and with how many columns.
/// @tparam Tout
/// @tparam Tin
/// @tparam Top
/// @param out
/// @param in
/// @param op
/// @return mods pointers
template <class Tout, class Tin, class Top>
QUICC_CUDA_HOSTDEV inline View::ViewBase<typename Tin::IndexType> getModsPointers(
   const Tout& out, const Tin& in, const Top& op)
{
   using namespace QuICC::View;
   using IndexType = typename Tin::IndexType;
   ViewBase<IndexType> modsPointers;
   using opLevelType = typename Top::LevelType;
   // JW uniform truncation projector/integrator ijk order
   if constexpr (std::is_same_v<opLevelType, CSL3D::level>)
   {
      modsPointers = in.pointers()[1];
   }
   // else if constexpr (ALegendre::is_integrator_v<Top>)
   else if constexpr (std::is_same_v<opLevelType, S1CLCSC3D::level> ||
                      std::is_same_v<opLevelType, S1CLCSC3DJIK::level>)
   {
      modsPointers = out.pointers()[1];
   }
   // else if constexpr (ALegendre::is_projector_v<Top>)
   else if constexpr (std::is_same_v<opLevelType, CS1RL3D::level> ||
                      std::is_same_v<opLevelType, CS1RL3DJIK::level>)
   {
      modsPointers = in.pointers()[1];
   }
   else
   {
      static_assert(std::is_same_v<opLevelType, void>,
         "backend for these types is not implemented.");
   }
   return modsPointers;
}


/// @brief POD to return matmul dimensions as a single object
/// @tparam T
template <class T> struct matmulDims
{
   T M, K, N;
};

/// @brief Get the matmul dimensions for the current layer
/// (col and layerIndex are layer dependent).
/// @tparam Tout
/// @tparam Tin
/// @tparam Top
/// @param out
/// @param in
/// @param op
/// @param col number of columns in the layer
/// @param layerIndex layer index
/// @return matmulDims
template <class Tout, class Tin, class Top>
QUICC_CUDA_HOSTDEV inline matmulDims<typename Tin::IndexType> getMatmulDims(
   const Tout& out, const Tin& in, const Top& op,
   const typename Tin::IndexType col, const typename Tin::IndexType layerIndex)
{
   using namespace QuICC::View;
   using IndexType = typename Tin::IndexType;
   using opLevelType = typename Top::LevelType;

   IndexType M, N, K;
   // JW uniform truncation projector/integrator
   if constexpr (std::is_same_v<opLevelType, CSL3D::level>)
   {
      M = op.dims()[0]; // radial points/modes
      K = op.dims()[1]; // radial modes/points
      N = col;          // number of columns
   }
   // else if constexpr (ALegendre::is_integrator_v<Top>)
   else if constexpr (std::is_same_v<opLevelType, S1CLCSC3D::level> ||
                      std::is_same_v<opLevelType, S1CLCSC3DJIK::level>)
   {
      IndexType L = out.dims()[0]; // harmonic order
      M = L - layerIndex;          // current harmonic order
      K = in.dims()[0];            // longitudinal points
      N = col;                     // number of columns
   }
   // else if constexpr (ALegendre::is_projector_v<Top>)
   else if constexpr (std::is_same_v<opLevelType, CS1RL3D::level> ||
                      std::is_same_v<opLevelType, CS1RL3DJIK::level>)
   {
      IndexType L = in.dims()[0]; // harmonic order
      M = out.dims()[0];          // longitudinal points
      K = L - layerIndex;         // current harmonic order
      N = col;                    // number of columns
   }
   else
   {
      static_assert(std::is_same_v<opLevelType, void>,
         "backend for these types is not implemented.");
   }

   return matmulDims<IndexType>{M, K, N};
}

} // namespace Quadrature
} // namespace Transform
} // namespace QuICC
