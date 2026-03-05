/**
 * @file SplittingAlgorithmDetails.hpp
 * @brief Base of the implementation of the load splitting algorithms
 */

#ifndef QUICC_PARALLEL_SPLITTINGALGORITHM_DETAILS_HPP
#define QUICC_PARALLEL_SPLITTINGALGORITHM_DETAILS_HPP

// System includes
//
#include <set>

// Project includes
//
#include "Profiler/Interface.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"

namespace QuICC {

namespace Parallel {

namespace details {

/// @brief Build 3D coordinates vector
template <Dimensions::Data::Id TDAT1D>
void getCoo(Dimensions::Transform::Id dimId, SharedResolution spRes,
   std::vector<std::array<int, 3>>& coos);

/// @brief Build 2D communication structure
void buildCommunicationStructure2D(const int localId, SharedResolution spRes,
   std::map<Dimensions::Transform::Id, std::multimap<int, int>>& commStructure);

/// @brief Build 3D communication structure
void buildCommunicationStructure3D(const int localId, SharedResolution spRes,
   std::map<Dimensions::Transform::Id, std::multimap<int, int>>& commStructure);

/// @brief Get global communication pattern
void getGlobalCommPattern(const int localId, const int nCpu,
   std::set<std::pair<int, int>>& filter);

template <Dimensions::Data::Id TDAT1D>
void getCoo(Dimensions::Transform::Id dimId, SharedResolution spRes,
   std::vector<std::array<int, 3>>& coos)
{
   Profiler::RegionFixture<4> fix(
      "Framework::LoadSplitter::SplittingAlgorithm::details::getCoo");

   const auto& tRes = *spRes->cpu()->dim(dimId);
   // Loop over third dimension
   for (int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
   {
      int k_ = tRes.idx<Dimensions::Data::DAT3D>(k);

      // Loop over second dimension
      for (int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
      {
         int j_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);

         // Loop over first dimension
         for (int i = 0; i < tRes.dim<TDAT1D>(j, k); i++)
         {
            int i_ = tRes.idx<TDAT1D>(i, j, k);

            // Generate point information
            auto point = spRes->counter().makeKey(dimId, i_, j_, k_);
            std::array<int, 3> p = {std::get<0>(point), std::get<1>(point),
               std::get<2>(point)};
            coos.push_back(p);
         }
      }
   }

   std::sort(coos.begin(), coos.end());
}

} // namespace details
} // namespace Parallel
} // namespace QuICC

#endif // QUICC_PARALLEL_SPLITTINGALGORITHM_DETAILS_HPP
