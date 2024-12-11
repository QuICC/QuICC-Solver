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


template <class TDATA>
void pack(View::ViewBase<TDATA> buffer, const TDATA* in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::Attributes<View::dense2D>> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView);



template <class TDATA>
void unPack(TDATA* out, const View::ViewBase<TDATA> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::Attributes<View::dense2D>> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView)

} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
