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
/// @param in
/// @param sendCountsView count for groupSize variables
/// @param sendDisplsView displacement for a single variable
/// @param sendBufferDisplsView displacement for groupSize variables
template <class TDATA>
void pack(View::ViewBase<TDATA> buffer, const std::vector<TDATA*> in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView);

/// @brief Unpack data from buffer, grouped version
/// @tparam TDATA
/// @param out
/// @param buffer
/// @param recvCountsView count for groupSize variables
/// @param recvDisplsView displacement for a single variable
/// @param recvBufferDisplsView displacement for groupSize variables
template <class TDATA>
void unPack(std::vector<TDATA*> out, const View::ViewBase<TDATA> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView);


} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
