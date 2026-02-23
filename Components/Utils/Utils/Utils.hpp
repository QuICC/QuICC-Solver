/**
 * @file Utils.hpp
 * @brief General utilities for (mainly) MPI enabled operations
 */
#pragma once

// External includes
//
#include <array>
#include <cassert>
#include <vector>

// Project includes
//
#include "Utils/Types.hpp"

namespace QuICC {
namespace Utils {

/// @brief Extract split index information from point_t indexes
void getSplitIdx(std::vector<int>& splitIdx, const std::vector<point_t>& coos);

/// @brief Match local and remote split index information
void matchSplitIdx(std::vector<int>& remNeededIdx,
   std::vector<int>& remNeededSizes, const std::vector<int>& locSplitIdx,
   const std::vector<int>& remSplitIdx);

/// @brief Filter point_t indexes based on split index information
void filterIdx(std::vector<point_t>& cooFiltered, std::vector<int>& cooSizes,
   const std::vector<point_t>& cooNew, const std::vector<int>& locIdx,
   const std::vector<int>& locDispl);

/// @brief Match local and remote point_t indexes to extract send displacements
void matchSendDispl(std::vector<std::vector<int>>& sendDispl,
   const std::vector<point_t>& locIdx, const std::vector<point_t>& remIdx,
   const std::vector<int>& remSizes, const std::vector<int>& remDispl);

/// @brief Match local and remote point_t indexes to extract send displacements
void matchSendDispl(std::vector<std::vector<int>>& sendDispl,
   const std::vector<point_t>& locIdx, const std::vector<point_t>& remIdx);

} // namespace Utils
} // namespace QuICC
