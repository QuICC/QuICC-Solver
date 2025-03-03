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
/// @brief namespace for Cpu backends
namespace Cpu {

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
   const View::ViewBase<int> sendBufferDisplsView)
{
   const auto I = sendDisplsView.dims()[0];
   const auto J = sendDisplsView.dims()[1];

   for (size_t i = 0; i < I; ++i)
   {
      for (int j = 0; j < sendCountsView[i]; ++j)
      {
         buffer[sendBufferDisplsView[i] + j] =
            *(in + sendDisplsView[i * J + j]);
      }
   }
}


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
   const View::ViewBase<int> recvBufferDisplsView)
{

   const auto I = recvDisplsView.dims()[0];
   const auto J = recvDisplsView.dims()[1];

   for (size_t i = 0; i < I; ++i)
   {
      for (int j = 0; j < recvCountsView[i]; ++j)
      {
         *(out + recvDisplsView[i * J + j]) =
            buffer[recvBufferDisplsView[i] + j];
      }
   }
}

} // namespace Cpu
} // namespace Transpose
} // namespace QuICC
