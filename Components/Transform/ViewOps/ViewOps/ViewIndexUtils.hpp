/**
 * @file ViewIndexUtils.hpp
 * @brief Utilities related data oredering and to building serial meta data
 */
#pragma once

// External includes
//
#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

// Project includes
//
#include "View/Attributes.hpp"
#include "View/ViewUtils.hpp"

namespace QuICC {
namespace View {
namespace Index {

/// @brief logical dimension ordering
/// these are used to identify that logical data ordering in a view
enum class dimOrder
{
   /// @brief xml file dimensions, logical ordering
   v012 = 0,
   /// @brief xml file dimensions
   RThetaPhi = 0,
   /// @brief implicit coriolis spectral spherical ordering
   NLM = 0,
   /// @brief xml file dimensions, cartesian ordering
   ZXY = 0,

   /// @brief logical ordering of dimensions
   v021 = 1,
   /// @brief physical space spherical ordering
   /// input to JW projector, output to AL integrator
   RPhiTheta = 1,
   /// @brief spectral spherical ordering
   /// input to JW projector, output to AL integrator
   NML = 1,
   /// @brief cartesian ordering of dimensions
   ZYX = 1,

   /// @brief logical ordering of dimensions
   v102 = 2,
   /// @brief physical space spherical ordering
   /// input to AL projector, output to FT integrator
   ThetaRPhi = 2,
   /// @brief spectral spherical ordering
   /// input to AL projector, output to FT integrator
   LNM = 2,
   /// @brief cartesian ordering of dimensions
   XZY = 2,

   /// @brief logical ordering of dimensions
   v210 = 3,
   /// @brief physical spherical ordering
   /// input to FT integrator, output to AL projector
   PhiThetaR = 3,
   /// @brief cartesian ordering of dimensions
   YXZ = 3,
};

/// @brief Reorder logical dimensions based on order
/// @param logDims logical dimensions, dimOrder::v012
/// @param order
/// @return
inline std::array<std::uint32_t, 3> getDims(
   const std::array<std::uint32_t, 3> logDims, const dimOrder order)
{
   std::array<std::uint32_t, 3> dims;
   switch (order)
   {
   case dimOrder::v012:
      dims = {logDims[0], logDims[1], logDims[2]};
      break;
   case dimOrder::v021:
      dims = {logDims[0], logDims[2], logDims[1]};
      break;
   case dimOrder::v102:
      dims = {logDims[1], logDims[0], logDims[2]};
      break;
   case dimOrder::v210:
      dims = {logDims[2], logDims[1], logDims[0]};
      break;
   default:
      throw std::logic_error("case not supported");
   }
   return dims;
}

/// @brief compute pointers and index meta data for fully populated
/// tensor at different stages based on view attributes
/// @tparam TagViewAttributes
/// @param logDims logical dimensions, dimOrder::v012
/// @param order order for which the metadata is created
/// @return pointer and indices in a struct
template <class TagViewAttributes>
ptrAndIdx densePtrAndIdx(const std::array<std::uint32_t, 3> logDims,
   const dimOrder order = dimOrder::v012)
{
   throw std::logic_error("not implemented");
   return ptrAndIdx{};
}


/// @brief compute pointers and index meta data for fully populated
/// tensor at different stages based on view attributes
/// @param logDims logical dimensions, dimOrder::v012
/// @param order order for which the metadata is created
/// @return pointer and indices in a struct
template <>
ptrAndIdx densePtrAndIdx<DCCSC3D>(const std::array<std::uint32_t, 3> logDims,
   const dimOrder order)
{
   std::array<std::uint32_t, 3> dims = getDims(logDims, order);
   auto N = dims[1];
   auto K = dims[2];
   // Populate meta for fully populated tensor
   ptrAndIdx ret;
   ret.ptr.resize(K + 1);
   ret.idx.resize(K * N);
   ret.ptr[0] = 0;
   for (std::size_t i = 1; i < ret.ptr.size(); ++i)
   {
      ret.ptr[i] = ret.ptr[i - 1] + N;
   }
   for (std::size_t i = 0; i < ret.idx.size(); ++i)
   {
      ret.idx[i] = i % N;
   }
   return ret;
}

/// @brief compute pointers and index meta data for fully populated
/// tensor at different stages based on view attributes
/// @param logDims logical dimensions, dimOrder::v012
/// @param order order for which the metadata is created
/// @return pointer and indices in a struct
template <>
ptrAndIdx densePtrAndIdx<DCCSC3DJIK>(const std::array<std::uint32_t, 3> logDims,
   const dimOrder order)
{
   return densePtrAndIdx<DCCSC3D>(logDims, order);
}

/// @brief compute pointers and index meta data for fully populated
/// tensor at different stages based on view attributes
/// @param logDims logical dimensions, dimOrder::v012
/// @param order order for which the metadata is created
/// @return pointer and indices in a struct
template <>
ptrAndIdx densePtrAndIdx<S1CLCSC3D>(const std::array<std::uint32_t, 3> logDims,
   const dimOrder order)
{
   return densePtrAndIdx<DCCSC3D>(logDims, order);
}

/// @brief compute pointers and index meta data for fully populated
/// tensor at different stages based on view attributes
/// @param logDims logical dimensions, dimOrder::v012
/// @param order order for which the metadata is created
/// @return pointer and indices in a struct
template <>
ptrAndIdx densePtrAndIdx<S1CLCSC3DJIK>(const std::array<std::uint32_t, 3> logDims,
   const dimOrder order)
{
   return densePtrAndIdx<S1CLCSC3D>(logDims, order);
}


/// @brief compute pointers and index meta data for fully populated
/// tensor at different stages based on view attributes
/// @tparam TagViewAttributes
/// @param logDims logical dimensions, dimOrder::v012
/// @param order order for which the metadata is created
/// @return pointer and indices in a struct
template <class TagViewAttributes>
ptrAndIdx densePtrAndIdxStep1(const std::array<std::uint32_t, 3> logDims,
   const dimOrder order = dimOrder::v012)
{
   throw std::logic_error("not implemented");
   return ptrAndIdx{};
}


/// @brief compute pointers and index meta data for fully populated
/// tensor at different stages based on view attributes
/// @param logDims logical dimensions, dimOrder::v012
/// @param order order for which the metadata is created
/// @return pointer and indices in a struct
template <>
ptrAndIdx densePtrAndIdxStep1<DCCSC3D>(
   const std::array<std::uint32_t, 3> logDims, const dimOrder order)
{
   std::array<std::uint32_t, 3> dims = getDims(logDims, order);
   auto K = dims[1];
   auto I = dims[2];

   // Populate meta for fully populated tensor
   ptrAndIdx ret;
   // row width (with jki) K - ...
   std::vector<std::uint32_t> kLoc(I);
   kLoc[0] = 1;
   for (std::size_t i = 1; i < I; ++i)
   {
      kLoc[i] = std::min(kLoc[i-1] + 1, K);
   }

   std::size_t layWidthCum = 0;
   ret.ptr.resize(I + 1);
   ret.ptr[0] = 0;
   for (std::size_t i = 1; i < I + 1; ++i)
   {
      std::uint32_t width = kLoc[i-1];
      ret.ptr[i] = ret.ptr[i - 1] + width;
      layWidthCum += width;
   }

   std::size_t layerIdx = 0;
   ret.idx.resize(layWidthCum);
   for (std::size_t l = 0; l < ret.ptr.size() - 1; ++l)
   {
      auto layerSize = ret.ptr[l + 1] - ret.ptr[l];
      for (std::size_t i = 0; i < layerSize; ++i)
      {
         ret.idx[layerIdx + i] = i;
      }
      layerIdx += layerSize;
   }
   return ret;
}

/// @brief compute pointers and index meta data for fully populated
/// tensor at different stages based on view attributes
/// @param logDims logical dimensions, dimOrder::v012
/// @param order order for which the metadata is created
/// @return pointer and indices in a struct
template <>
ptrAndIdx densePtrAndIdxStep1<DCCSC3DJIK>(const std::array<std::uint32_t, 3> logDims,
   const dimOrder order)
{
   return densePtrAndIdxStep1<DCCSC3D>(logDims, order);
}

} // namespace Index
} // namespace View
} // namespace QuICC
