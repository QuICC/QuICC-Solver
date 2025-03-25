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

   for (std::size_t i = 0; i < I; ++i)
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

   for (std::size_t i = 0; i < I; ++i)
   {
      for (int j = 0; j < recvCountsView[i]; ++j)
      {
         *(out + recvDisplsView[i * J + j]) =
            buffer[recvBufferDisplsView[i] + j];
      }
   }
}

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
   const View::ViewBase<int> sendBufferDisplsView, const std::int64_t groupSize)
{
   const int I = sendDisplsView.dims()[0];
   const int J = sendDisplsView.dims()[1];

   for (int g = 0; g < groupSize; ++g)
   {
      for (int i = 0; i < I; ++i)
      {
         int sendCount = sendCountsView[i] / groupSize;
         for (int j = 0; j < sendCount; ++j)
         {
            buffer[g * sendCount + sendBufferDisplsView[i] + j] =
               *(in[g] + sendDisplsView[i * J + j]);
         }
      }
   }
}

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
   const View::ViewBase<int> recvBufferDisplsView, const std::int64_t groupSize)
{

   const int I = recvDisplsView.dims()[0];
   const int J = recvDisplsView.dims()[1];

   for (int g = 0; g < groupSize; ++g)
   {
      for (int i = 0; i < I; ++i)
      {
         int recvCount = recvCountsView[i] / groupSize;
         for (int j = 0; j < recvCount; ++j)
         {
            *(out[g] + recvDisplsView[i * J + j]) =
               buffer[g * recvCount + recvBufferDisplsView[i] + j];
         }
      }
   }
}

} // namespace Cpu
} // namespace Transpose
} // namespace QuICC
