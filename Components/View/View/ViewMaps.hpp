/**
 * @file ViewMaps.hpp
 * @brief Affine mappings to access dense indices
 */

#pragma once

// System includes
//
#include <cassert>
#include <stdexcept>
#include <vector>

// Project includes
//
#include "View/ViewMacros.hpp"



namespace QuICC {
namespace Memory {

/// @brief Map in-order index (complex DFT) [-N...N] -> [0...N,-N...-1]
/// with padding
/// @tparam IndexType
/// @param i logical index to be mappe
/// @param N logical dense dimension
/// @param NPadded in memory dense dimension
/// @return in memory index
template <class IndexType>
inline QUICC_CUDA_HOST IndexType mapInOrderIndex(const IndexType i,
   const IndexType N, const IndexType NPadded)
{
   const auto shift = std::round(static_cast<double>(i)/N)* ((N+1)/2 + NPadded-N);
   const auto base = i % ((N+1) / 2);
   return base + shift;
}

/// @brief Map in-order index (complex DFT) [-N...N] -> [0...N,-N...-1]
/// without padding
/// @tparam IndexType
/// @param i logical index to be mappe
/// @param N logical dense dimension
/// @return in memory index
template <class IndexType>
inline QUICC_CUDA_HOST IndexType mapInOrderIndex(const IndexType i,
   const IndexType N)
{
   return mapInOrderIndex(i, N, N);
}

} // namespace Memory
} // namespace QuICC
