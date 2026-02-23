/**
 * @file Utils.hpp
 * @brief General utilities for (mainly) MPI enabled operations
 */
#pragma once

// External includes
//
#include <array>
#include <cassert>
#include <memory>
#include <mpi.h>
#include <vector>

// Project includes
//
#include "Utils/Types.hpp"

namespace QuICC {
namespace Utils {
namespace Mpi {

/// @brief Broadcast initial split index information
const std::vector<int>* broadcastSplitIdx(const int r, const int rank,
   std::vector<int>& remSplit, const std::vector<int>& locSplit,
   const MPI_Comm comm);

/// @brief Distribute sizes of split indexes to all ranks
void distributeSplitSizes(const std::vector<int>& sendSizes,
   std::vector<int>& sendDispl, std::vector<int>& recvSizes,
   std::vector<int>& recvDispl, const MPI_Comm comm);

/// @brief Distribute split indexes to all ranks
void distributeSplitIdx(const std::vector<int>& sendIdx,
   const std::vector<int>& sendSizes, std::vector<int>& recvIdx,
   std::vector<int>& recvSizes, std::vector<int>& recvDispl,
   const MPI_Comm comm);

/// @brief Distribute split point_t indexes to all ranks
void distributeSplitIdx(const std::vector<point_t>& sendIdx,
   const std::vector<int>& sendSizes, std::vector<point_t>& recvIdx,
   std::vector<int>& recvSizes, std::vector<int>& recvDispl,
   const MPI_Comm comm);

} // namespace Mpi
} // namespace Utils
} // namespace QuICC
