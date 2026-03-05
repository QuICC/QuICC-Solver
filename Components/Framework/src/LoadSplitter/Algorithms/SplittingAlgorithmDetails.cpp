/**
 * @file SplittingAlgorithmDetails.cpp
 * @brief Source of the base of the implementation of the load splitting
 * algorithms
 */

// System includes
//
#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>

#ifdef QUICC_MPI
#include <mpi.h>
#endif

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Enums/DimensionTools.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithmDetails.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "Utils/Utils.hpp"
#ifdef QUICC_MPI
#include "Utils/Mpi/Utils.hpp"
#endif

namespace QuICC {

namespace Parallel {

namespace details {

void buildCommunicationStructure2D(const int localId, SharedResolution spRes,
   std::map<Dimensions::Transform::Id, std::multimap<int, int>>& commStructure)
{
   Profiler::RegionFixture<4> fix("Framework::LoadSplitter::SplittingAlgorithm:"
                                  ":details::buildCommunicationStructure2D");

   Dimensions::Transform::Id dimId;
   int i_;
   int j_;
   int k_;

   // Simplify syntax
   typedef std::pair<int, int> Coordinate;

   // Extract communication structure from resolution object
   std::set<Coordinate> bwdMap;
   std::set<Coordinate> fwdMap;

   // Storage for a coordinate
   Coordinate point;

   // Position iterator for insert calls
   std::set<Coordinate>::iterator mapPos;

   // Loop over possible data exchanges
   std::vector<Dimensions::Transform::Id> exchanges = {
      Dimensions::Transform::TRA1D};
   for (auto exId: exchanges)
   {
      // Create storage for structure
      commStructure.emplace(exId, std::multimap<int, int>());

      // initialise the position hint for inserts
      mapPos = bwdMap.begin();

      dimId = Dimensions::jump(exId, 1);
      const auto& bwdTRes = *spRes->cpu()->dim(dimId);
      // Loop over second dimension
      for (int j = 0; j < bwdTRes.dim<Dimensions::Data::DAT2D>(); j++)
      {
         j_ = bwdTRes.idx<Dimensions::Data::DAT2D>(j);

         // Loop over backward dimension
         for (int i = 0; i < bwdTRes.dim<Dimensions::Data::DATB1D>(); i++)
         {
            i_ = bwdTRes.idx<Dimensions::Data::DATB1D>(i);

            // Generate point information
            point = spRes->counter().makeKey(dimId, i_, j_);

            // Get insertion position to use as next starting point to speed up
            // insertion
            mapPos = bwdMap.insert(mapPos, point);
         }
      }

      // Loop over CPUs
      MatrixI matRemote;
      int matched = 0;
      int toMatch = -1;
      std::set<std::pair<int, int>> filter;
      dimId = exId;
      const auto& fwdTRes = *spRes->cpu()->dim(dimId);
      for (int cpu = 0; cpu < spRes->nCpu(); cpu++)
      {
         matched = 0;

         // Local CPU
         if (cpu == localId)
         {
            // Loop over second dimension
            for (int j = 0; j < fwdTRes.dim<Dimensions::Data::DAT2D>(); j++)
            {
               j_ = fwdTRes.idx<Dimensions::Data::DAT2D>(j);

               // Loop over forward dimension
               for (int i = 0; i < fwdTRes.dim<Dimensions::Data::DATF1D>(); i++)
               {
                  i_ = fwdTRes.idx<Dimensions::Data::DATF1D>(i);

                  // Generate point information
                  point = spRes->counter().makeKey(dimId, i_, j_);

                  // Look for same key in backward list
                  mapPos = bwdMap.find(point);

                  // Key was present, drop entry and extend filter
                  if (mapPos != bwdMap.end())
                  {
                     // Add corresponding communication edge to filter
                     filter.insert(std::make_pair(cpu, localId));

                     // Delete found coordinate
                     bwdMap.erase(mapPos);
                  }
                  else
                  {
                     fwdMap.insert(point);
                  }
               }
            }

            // Store size of forward coordinates
            toMatch = fwdMap.size();

#ifdef QUICC_MPI
            // Convert coordinates set to matrix to send through MPI
            matRemote.resize(2, fwdMap.size());
            int i = 0;
            for (auto it = fwdMap.begin(); it != fwdMap.end(); ++it)
            {
               matRemote(0, i) = it->first;
               matRemote(1, i) = it->second;
               i++;
            }

            // Broadcast size
            QuICCEnv().synchronize();
            int ierr = MPI_Bcast(&toMatch, 1, MPI_INT, cpu, MPI_COMM_WORLD);
            QuICCEnv().check(ierr, 711);

            // Broadcast data
            QuICCEnv().synchronize();
            ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpu,
               MPI_COMM_WORLD);
            QuICCEnv().check(ierr, 712);

            // Remote CPU
         }
         else
         {
            // Get size
            QuICCEnv().synchronize();
            int ierr = MPI_Bcast(&toMatch, 1, MPI_INT, cpu, MPI_COMM_WORLD);
            QuICCEnv().check(ierr, 713);

            // Get remote keys as matrix
            matRemote.resize(2, toMatch);
            QuICCEnv().synchronize();
            ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpu,
               MPI_COMM_WORLD);
            QuICCEnv().check(ierr, 714);

            // Compare received data to stored indexes
            for (int i = 0; i < toMatch; i++)
            {
               point = std::make_pair(matRemote(0, i), matRemote(1, i));

               mapPos = bwdMap.find(point);

               // Check if point is in backward map
               if (mapPos != bwdMap.end())
               {
                  // Add corresponding communication edge to filter
                  filter.insert(std::make_pair(cpu, localId));

                  // Delete found entry
                  bwdMap.erase(mapPos);

                  // Increase matched counter
                  matched++;
               }
            }
         }

#else
         }
#endif // QUICC_MPI
      }

#ifdef QUICC_MPI
      // Gather total number of match entries
      QuICCEnv().synchronize();
      int ierr = MPI_Allreduce(MPI_IN_PLACE, &matched, 1, MPI_INT, MPI_SUM,
         MPI_COMM_WORLD);
      QuICCEnv().check(ierr, 715);
#endif // QUICC_MPI

      // Check that everything matched
      if (toMatch != matched)
      {
         throw std::logic_error("The computed index sets don't match!");
      }

      // Get global minimized communcation pattern
      details::getGlobalCommPattern(localId, spRes->nCpu(), filter);

      // Store obtained minimized structure
      for (auto&& fId: filter)
      {
         commStructure.at(exId).insert(fId);
      }

      // Clear all the data for next loop
      bwdMap.clear();
      fwdMap.clear();
   }

   // Synchronize
   QuICCEnv().synchronize();
}

void buildCommunicationStructure3D(const int localId, SharedResolution spRes,
   std::map<Dimensions::Transform::Id, std::multimap<int, int>>& commStructure)
{
   Profiler::RegionFixture<4> fix("Framework::LoadSplitter::SplittingAlgorithm:"
                                  ":details::buildCommunicationStructure3D");

   Dimensions::Transform::Id dimId;
   int i_;
   int j_;
   int k_;

   // Simplify syntax
   typedef std::array<int, 3> point_t;

   // Loop over possible data exchanges
   std::vector<Dimensions::Transform::Id> exchanges = {
      Dimensions::Transform::TRA1D, Dimensions::Transform::TRA2D,
      Dimensions::Transform::SPECTRAL};
   for (auto exId: exchanges)
   {
      // Create storage for structure
      commStructure.emplace(exId, std::multimap<int, int>());

      dimId = Dimensions::jump(exId, 1);

      std::vector<point_t> absCooOld;
      std::vector<point_t> absCooNew;
      details::getCoo<Dimensions::Data::DATB1D>(dimId, spRes, absCooNew);

      details::getCoo<Dimensions::Data::DATF1D>(exId, spRes, absCooOld);

#ifdef QUICC_MPI
      int ranks = spRes->nCpu();
      int rank = localId;
      MPI_Comm comm = MPI_COMM_WORLD;

      std::vector<int> locOldIdxSplit;
      Utils::getSplitIdx(locOldIdxSplit, absCooOld);

      std::vector<int> remNewIdxSplit;
      Utils::getSplitIdx(remNewIdxSplit, absCooNew);

      std::vector<int> remNewIdxSplitNeededAll;
      std::vector<int> remNewIdxSplitNeededSizes;
      for (int r = 0; r < ranks; ++r)
      {
         // get new coo from other rank and check if it is here
         std::vector<int> tmpSplit;
         const std::vector<int>* pRemSplit;

         pRemSplit = Utils::Mpi::broadcastSplitIdx(r, rank, tmpSplit,
            remNewIdxSplit, comm);

         // Match split indexes
         Utils::matchSplitIdx(remNewIdxSplitNeededAll,
            remNewIdxSplitNeededSizes, locOldIdxSplit, *pRemSplit);
      }

      // Distributed split indexes
      std::vector<int> locNewIdxSplitNeededAll;
      std::vector<int> locNewIdxSplitNeededSizes;
      std::vector<int> locNewIdxSplitNeededDispl;
      Utils::Mpi::distributeSplitIdx(remNewIdxSplitNeededAll,
         remNewIdxSplitNeededSizes, locNewIdxSplitNeededAll,
         locNewIdxSplitNeededSizes, locNewIdxSplitNeededDispl, comm);

      // Filter new indexes
      std::vector<point_t> absCooNewAll;
      std::vector<int> absCooNewSizes;
      Utils::filterIdx(absCooNewAll, absCooNewSizes, absCooNew,
         locNewIdxSplitNeededAll, locNewIdxSplitNeededDispl);

      // Distributed filtered split indexes
      std::vector<point_t> remAbsCooNewAll;
      std::vector<int> remAbsCooNewSizes;
      std::vector<int> remAbsCooNewDispl;
      Utils::Mpi::distributeSplitIdx(absCooNewAll, absCooNewSizes,
         remAbsCooNewAll, remAbsCooNewSizes, remAbsCooNewDispl, comm);

      // Find matches for sendDispls
      std::vector<std::vector<int>> sendDispls;
      Utils::matchSendDispl(sendDispls, absCooOld, remAbsCooNewAll,
         remAbsCooNewSizes, remAbsCooNewDispl);
#else
      std::vector<std::vector<int>> sendDispls;
      Utils::matchSendDispl(sendDispls, absCooOld, absCooNew);
#endif // QUICC_MPI

      int toMatch = absCooOld.size();
      int matched = 0;
      std::set<std::pair<int, int>> filter;
      for (int r = 0; r < sendDispls.size(); r++)
      {
         if (sendDispls.at(r).size() > 0)
         {
            filter.emplace(r, localId);
            matched += sendDispls.at(r).size();
         }
      }

      // Check that everything matched
      if (toMatch != matched)
      {
         throw std::logic_error("The computed index sets don't match!");
      }

      // Get global minimized communcation pattern
      details::getGlobalCommPattern(localId, spRes->nCpu(), filter);

      // Store obtained minimized structure
      for (auto&& fId: filter)
      {
         commStructure.at(exId).insert(fId);
      }
   }

   // Synchronize
   QuICCEnv().synchronize();
}

void getGlobalCommPattern(const int localId, const int nCpu,
   std::set<std::pair<int, int>>& filter)
{
#ifdef QUICC_MPI
   Profiler::RegionFixture<4> fix("Framework::LoadSplitter::SplittingAlgorithm:"
                                  ":details::getGlobalCommPattern");

   std::vector<std::array<int, 2>> commFilter;

   // Gather full communication structure
   for (int cpu = 0; cpu < nCpu; cpu++)
   {
      int filterSize = 0;
      if (cpu == localId)
      {
         // Send local filter
         commFilter.reserve(filter.size());
         for (auto it = filter.begin(); it != filter.end(); ++it)
         {
            std::array<int, 2> c = {it->first, it->second};
            commFilter.emplace_back(c);
         }

         filterSize = 2 * commFilter.size();

         // Get size
         QuICCEnv().synchronize();
         int ierr = MPI_Bcast(&filterSize, 1, MPI_INT, cpu, MPI_COMM_WORLD);
         QuICCEnv().check(ierr, 725);

         // Get remote comm pattern
         QuICCEnv().synchronize();
         ierr = MPI_Bcast(commFilter.data(), 2 * commFilter.size(), MPI_INT,
            cpu, MPI_COMM_WORLD);
         QuICCEnv().check(ierr, 726);
      }
      else
      {
         // Get size
         QuICCEnv().synchronize();
         int ierr = MPI_Bcast(&filterSize, 1, MPI_INT, cpu, MPI_COMM_WORLD);
         QuICCEnv().check(ierr, 727);

         // Get remote comm pattern
         commFilter.resize(filterSize / 2);
         QuICCEnv().synchronize();
         ierr = MPI_Bcast(commFilter.data(), filterSize, MPI_INT, cpu,
            MPI_COMM_WORLD);
         QuICCEnv().check(ierr, 728);

         for (auto it = commFilter.begin(); it != commFilter.end(); ++it)
         {
            filter.insert(std::make_pair((*it)[0], (*it)[1]));
         }
      }
   }
#endif // QUICC_MPI
}

} // namespace details
} // namespace Parallel
} // namespace QuICC
