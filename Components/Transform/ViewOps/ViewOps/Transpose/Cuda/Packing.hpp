/**
 * @file Op.hpp
 * @brief Transpose operations on Views
 */
#pragma once

// External includes
//

// Project includes
//
#include "View/View.hpp"
#include "ViewOps/Transpose/StructArray.hpp"

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cuda backends
namespace Cuda {

/// @brief Pack data into buffer
/// @tparam TDATA
/// @param buffer
/// @param in
/// @param sendCountsView
/// @param sendDisplsView
/// @param sendBufferDisplsView
template <class TDATA>
void pack(View::ViewBase<TDATA> buffer, const TDATA* in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView);


/// @brief Unpack data from buffer
/// @tparam TDATA
/// @param out
/// @param buffer
/// @param recvCountsView
/// @param recvDisplsView
/// @param recvBufferDisplsView
template <class TDATA>
void unPack(TDATA* out, const View::ViewBase<TDATA> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView);


/// @brief Pack data into buffer, grouped version
/// @tparam TDATA
/// @param buffer
/// @param in array of pointers to data (size is fixed and >= groupSize)
/// @param sendCountsView count for groupSize variables
/// @param sendDisplsView displacement for a single variable
/// @param sendBufferDisplsView displacement for groupSize variables
/// @param groupSize actual group size
template <class TDATA, int SIZE>
void pack(View::ViewBase<TDATA> buffer, structArray<const TDATA*, SIZE> in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView, const std::int64_t groupSize);

/// @brief Unpack data from buffer, grouped version
/// @tparam TDATA
/// @param out array of pointers to data (size is fixed and >= groupSize)
/// @param buffer
/// @param recvCountsView count for groupSize variables
/// @param recvDisplsView displacement for a single variable
/// @param recvBufferDisplsView displacement for groupSize variables
/// @param groupSize actual group size
template <class TDATA, int SIZE>
void unPack(structArray<TDATA*, SIZE> out, const View::ViewBase<TDATA> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView, const std::int64_t groupSize);


} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
