/**
 * @file Utils.cpp
 * @brief General utilities
 */

// External includes
//
#include <algorithm>
#include <numeric>
#include <set>

// Project includes
//
#include "Profiler/Interface.hpp"
#include "Utils/Utils.hpp"

namespace QuICC {
namespace Utils {

void getSplitIdx(std::vector<int>& splitIdx, const std::vector<point_t>& coos)
{
   Profiler::RegionFixture<4> fix("Utils::getSplitIdx");

   const auto dimSize = std::tuple_size<point_t>{};

   // Extract index set for each dimension
   std::vector<std::set<int>> filters(dimSize);
   for (std::size_t i = 0; i < coos.size(); ++i)
   {
      auto&& p = coos[i];
      for (int i = 0; i < filters.size(); i++)
      {
         filters.at(i).insert(p[i]);
      }
   }

   // Compute total size
   int tot = 0;
   for (auto&& s: filters)
   {
      tot += s.size();
   }

   // Set pointers for accessing indexes for each dimensions
   splitIdx.reserve(filters.size() + dimSize + 1 + tot);
   splitIdx.push_back(filters.size() + 1);
   for (auto&& s: filters)
   {
      splitIdx.push_back(splitIdx.back() + s.size());
   }

   // Add all indexes
   for (auto&& s: filters)
   {
      std::copy(s.begin(), s.end(), std::back_inserter(splitIdx));
   }
   assert(splitIdx.size() == splitIdx.at(3));
}

void matchSplitIdx(std::vector<int>& remNeededIdx,
   std::vector<int>& remNeededSizes, const std::vector<int>& locSplitIdx,
   const std::vector<int>& remSplitIdx)
{
   Profiler::RegionFixture<4> fix("Utils::matchSplitIdx");

   const auto dimSize = std::tuple_size<point_t>{};
   // loop over loc coo to find match
   std::vector<int> remSplitIdxNeeded(dimSize + 1, 0);
   for (int i = 0; i < dimSize; i++)
   {
      remSplitIdxNeeded.at(i) = remSplitIdxNeeded.size();
      if (locSplitIdx.size() > dimSize + 1 && remSplitIdx.size() > dimSize + 1)
      {
         std::set_intersection(locSplitIdx.begin() + locSplitIdx.at(i),
            locSplitIdx.begin() + locSplitIdx.at(i + 1),
            remSplitIdx.begin() + remSplitIdx.at(i),
            remSplitIdx.begin() + remSplitIdx.at(i + 1),
            std::back_inserter(remSplitIdxNeeded));
      }
   }
   remSplitIdxNeeded.at(dimSize) = remSplitIdxNeeded.size();
   std::copy(remSplitIdxNeeded.begin(), remSplitIdxNeeded.end(),
      std::back_inserter(remNeededIdx));
   remNeededSizes.push_back(remSplitIdxNeeded.size());
}

void filterIdx(std::vector<point_t>& cooFiltered, std::vector<int>& cooSizes,
   const std::vector<point_t>& cooNew, const std::vector<int>& locIdx,
   const std::vector<int>& locDispl)
{
   Profiler::RegionFixture<4> fix("Utils::filterIdx");

   cooFiltered.reserve(cooNew.size());
   for (int r = 0; r < locDispl.size(); ++r)
   {
      int istart = cooFiltered.size();
      int count = 0;
      auto itStart = locIdx.begin() + locDispl.at(r);
      auto itStart0 = itStart + (*itStart);
      auto itEnd0 = itStart + (*(itStart + 1));
      auto itStart1 = itStart + (*(itStart + 1));
      auto itEnd1 = itStart + (*(itStart + 2));
      auto itStart2 = itStart + (*(itStart + 2));
      auto itEnd2 = itStart + (*(itStart + 3));
      for (auto&& p: cooNew)
      {
         if (std::binary_search(itStart0, itEnd0, p[0]))
         {
            if (std::binary_search(itStart1, itEnd1, p[1]))
            {
               if (std::binary_search(itStart2, itEnd2, p[2]))
               {
                  cooFiltered.push_back(p);
                  count++;
               }
            }
         }
      }
      cooSizes.push_back(count);
      std::sort(cooFiltered.begin() + istart, cooFiltered.end());
   }
}

void matchSendDispl(std::vector<std::vector<int>>& sendDispl,
   const std::vector<point_t>& locIdx, const std::vector<point_t>& remIdx,
   const std::vector<int>& remSizes, const std::vector<int>& remDispl)
{
   Profiler::RegionFixture<4> fix("Utils::matchSendDispl");

   // Compute permutation indexes
   std::vector<int> argsort(locIdx.size());
   std::iota(argsort.begin(), argsort.end(), 0);
   std::sort(argsort.begin(), argsort.end(),
      [&](std::size_t i, std::size_t j) { return locIdx[i] < locIdx[j]; });

   // Reserve memory
   sendDispl.resize(remDispl.size());
   for (int r = 0; r < remDispl.size(); ++r)
   {
      sendDispl[r].reserve(
         std::min(locIdx.size(), static_cast<std::size_t>(remSizes.at(r))));
   }

   // Extract matching list
   for (int i: argsort)
   {
      auto&& p = locIdx[i];

      for (int r = 0; r < remDispl.size(); ++r)
      {
         // get new coo from other rank and check if it is here
         auto remBegin = remIdx.begin() + remDispl.at(r);
         auto remEnd = remIdx.begin() + remDispl.at(r) + remSizes.at(r);
         assert(std::is_sorted(remBegin, remEnd));

         if (std::binary_search(remBegin, remEnd, p))
         {
            sendDispl[r].emplace_back(i);
            break;
         }
      }
   }
}

void matchSendDispl(std::vector<std::vector<int>>& sendDispl,
   const std::vector<point_t>& locIdx, const std::vector<point_t>& remIdx)
{
   Profiler::RegionFixture<4> fix("Utils::matchSendDisplSerial");

   // Compute permutation indexes
   std::vector<int> argsort(locIdx.size());
   std::iota(argsort.begin(), argsort.end(), 0);
   std::sort(argsort.begin(), argsort.end(),
      [&](std::size_t i, std::size_t j) { return locIdx[i] < locIdx[j]; });

   // Reserve memory
   sendDispl.resize(1);
   sendDispl[0].reserve(std::min(locIdx.size(), remIdx.size()));

   // Extract matching list
   for (int i: argsort)
   {
      auto&& p = locIdx[i];

      // get new coo from other rank and check if it is here
      assert(std::is_sorted(remIdx.begin(), remIdx.end()));

      if (std::binary_search(remIdx.begin(), remIdx.end(), p))
      {
         sendDispl[0].emplace_back(i);
      }
   }
}

} // namespace Utils
} // namespace QuICC
