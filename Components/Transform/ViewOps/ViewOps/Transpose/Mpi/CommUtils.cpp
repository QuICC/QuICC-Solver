/**
 * @file Comm.cpp
 * @brief Methods for mpi enabled transform
 */

// External includes
//
#include <algorithm>
#include <map>
#include <set>

// Project includes
//
#include "ViewOps/Transpose/Mpi/CommUtils.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Transpose {
namespace Mpi {

#if 0
std::vector<std::vector<int>> getDispls(const std::vector<point_t>& absCooNew,
   const std::vector<point_t>& absCooOld, const MPI_Comm comm)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::getDispls");
   int rank, ranks;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &ranks);

   std::vector<std::vector<int>> sendDispls(ranks);
   std::map<point_t, int> locOldIdx;
   for (std::size_t i = 0; i < absCooOld.size(); ++i)
   {
      auto&& p = absCooOld[i];
      locOldIdx[p] = i;
   }
   for (int r = 0; r < ranks; ++r)
   {
      // get new coo from other rank and check if it is here
      std::map<point_t, int> remNewIdx;
      // comm remote coo size
      int remAbsCooNewSize = absCooNew.size();
      MPI_Bcast(&remAbsCooNewSize, 1, MPI_INT, r, comm);
      if (r == rank)
      {
         MPI_Bcast(const_cast<point_t*>(absCooNew.data()),
            absCooNew.size() * dimSize, MPI_INT, r, comm);
         // setup remote map
         for (std::size_t i = 0; i < absCooNew.size(); ++i)
         {
            auto&& p = absCooNew[i];
            remNewIdx[p] = i;
         }
      }
      else
      {
         // comm remote coordinates
         std::vector<point_t> remAbsCooNew(remAbsCooNewSize);
         MPI_Bcast(remAbsCooNew.data(), remAbsCooNew.size() * dimSize, MPI_INT,
            r, comm);
         // setup remote map
         for (std::size_t i = 0; i < remAbsCooNew.size(); ++i)
         {
            auto&& p = remAbsCooNew[i];
            remNewIdx[p] = i;
         }
      }

      // loop over loc coo to find match
      for (auto itLCoo = locOldIdx.begin(); itLCoo != locOldIdx.end();)
      {
         auto lCoo = (*itLCoo).first;
         if (auto itRCoo = remNewIdx.find(lCoo); itRCoo != remNewIdx.end())
         {
            sendDispls[r].push_back((*itLCoo).second);
            itLCoo = locOldIdx.erase(itLCoo);
            remNewIdx.erase(itRCoo);
         }
         else
         {
            ++itLCoo;
         }
      }
   }
   return sendDispls;
}
#else

namespace details {
void getSplitIdx(std::vector<int>& splitIdx, const std::vector<point_t>& coos)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::details::getSplitIdx");

   // Extract index set for each dimension
   std::vector<std::set<int>> filters(std::tuple_size<point_t>{});
   for (std::size_t i = 0; i < coos.size(); ++i)
   {
      auto&& p = coos[i];
      for(int i = 0; i < filters.size(); i++)
      {
         filters.at(i).insert(p[i]);
      }
   }

   // Compute total size
   int tot = 0;
   for(auto&& s: filters)
   {
      tot += s.size();
   }

   // Set pointers for accessing indexes for each dimensions
   splitIdx.reserve(filters.size() + dimSize + 1 + tot);
   splitIdx.push_back(filters.size() + 1);
   for(auto&& s: filters)
   {
      splitIdx.push_back(splitIdx.back() + s.size());
   }

   // Add all indexes
   for(auto&& s: filters)
   {
      std::copy(s.begin(), s.end(), std::back_inserter(splitIdx));
   }
   assert(splitIdx.size() == splitIdx.at(3));
}

void matchSplitIdx(std::vector<int>& remNeededIdx, std::vector<int>& remNeededSizes, const std::vector<int>& locSplitIdx, const std::vector<int>& remSplitIdx)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::details::matchSplitIdx");

   // loop over loc coo to find match
   std::vector<int>  remSplitIdxNeeded(dimSize + 1, 0);
   for(int i = 0; i < dimSize; i++)
   {
      remSplitIdxNeeded.at(i) = remSplitIdxNeeded.size();
      if(locSplitIdx.size() > dimSize + 1 && remSplitIdx.size() > dimSize + 1)
      {
         std::set_intersection(
               locSplitIdx.begin() + locSplitIdx.at(i), locSplitIdx.begin() + locSplitIdx.at(i+1),
               remSplitIdx.begin() + remSplitIdx.at(i), remSplitIdx.begin() + remSplitIdx.at(i+1),
               std::back_inserter(remSplitIdxNeeded));
      }
   }
   remSplitIdxNeeded.at(dimSize) = remSplitIdxNeeded.size();
   std::copy(remSplitIdxNeeded.begin(), remSplitIdxNeeded.end(), std::back_inserter(remNeededIdx));
   remNeededSizes.push_back(remSplitIdxNeeded.size());
}

const std::vector<int>* broadcastSplitIdx(const int r, const int rank, std::vector<int>& remSplit, const std::vector<int>& locSplit, const MPI_Comm comm)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::details::broadcastSplitIdx");

   // get new coo from other rank and check if it is here
   const std::vector<int>* pRemSplit;

   if (r == rank)
   {
      pRemSplit = &locSplit;
   }
   else
   {
      pRemSplit = &remSplit;
   }

   // comm remote coo size
   int remAbsCooNewSize = pRemSplit->size();
   MPI_Bcast(&remAbsCooNewSize, 1, MPI_INT, r, comm);

   if (r == rank)
   {
      MPI_Bcast(const_cast<int*>(locSplit.data()),
            locSplit.size(), MPI_INT, r, comm);
   }
   else
   {
      // comm remote coordinates
      remSplit.resize(remAbsCooNewSize);
      MPI_Bcast(remSplit.data(), remSplit.size(), MPI_INT,
            r, comm);
   }

   return pRemSplit;
}

void distributeSplitSizes(const std::vector<int>& sendSizes, std::vector<int>& sendDispl, std::vector<int>& recvSizes, std::vector<int>& recvDispl, const MPI_Comm comm)
{
   Profiler::RegionFixture<5> fix("Transpose::Mpi::OpGrouped::details::distributeSplitSizes");

   int ranks;
   MPI_Comm_size(comm, &ranks);

   recvSizes.resize(ranks);
   MPI_Alltoall(sendSizes.data(), 1, MPI_INT, recvSizes.data(), 1, MPI_INT, comm);

   recvDispl.clear();
   sendDispl.clear();
   recvDispl.push_back(0);
   sendDispl.push_back(0);
   for(int i = 0; i < ranks-1; i++)
   {
      recvDispl.push_back(recvDispl.back() + recvSizes.at(i));
      sendDispl.push_back(sendDispl.back() + sendSizes.at(i));
   }
}

void distributeSplitIdx(const std::vector<int>& sendIdx, const std::vector<int>& sendSizes, std::vector<int>& recvIdx, std::vector<int>& recvSizes, std::vector<int>& recvDispl, const MPI_Comm comm)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::details::distributeSplitIdx");

   std::vector<int> sendDispl;
   distributeSplitSizes(sendSizes, sendDispl, recvSizes, recvDispl, comm);

   int tot = 0;
   for(auto&& s: recvSizes)
   {
      tot += s;
   }

   recvIdx.resize(tot);
   MPI_Alltoallv(sendIdx.data(), sendSizes.data(), sendDispl.data(), MPI_INT, recvIdx.data(), recvSizes.data(), recvDispl.data(), MPI_INT, comm);
}

void distributeSplitIdx(const std::vector<point_t>& sendIdx, const std::vector<int>& sendSizes, std::vector<point_t>& recvIdx, std::vector<int>& recvSizes, std::vector<int>& recvDispl, const MPI_Comm comm)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::details::distributeSplitIdx_point");

   std::vector<int> sendDispl;
   distributeSplitSizes(sendSizes, sendDispl, recvSizes, recvDispl, comm);

   int tot = 0;
   for(auto&& s: recvSizes)
   {
      tot += s;
   }
   recvIdx.resize(tot);

   std::vector<int> sendSizes_(sendSizes);
   std::vector<int> sendDispl_(sendDispl);
   std::vector<int> recvSizes_(recvSizes);
   std::vector<int> recvDispl_(recvDispl);
   for(int i = 0; i < sendDispl.size(); i++)
   {
      sendSizes_.at(i) *= dimSize;
      sendDispl_.at(i) *= dimSize;
      recvSizes_.at(i) *= dimSize;
      recvDispl_.at(i) *= dimSize;
   }

   MPI_Alltoallv(sendIdx.data(), sendSizes_.data(), sendDispl_.data(), MPI_INT, recvIdx.data(), recvSizes_.data(), recvDispl_.data(), MPI_INT, comm);
}

void filterIdx(std::vector<point_t>& cooFiltered, std::vector<int>& cooSizes, const std::vector<point_t>& cooNew, const std::vector<int>& locIdx, const std::vector<int>& locDispl)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::details::filterIdx");

   for(int r = 0; r < locDispl.size(); ++r)
   {
      int count = 0;
      auto itStart = locIdx.begin() + locDispl.at(r);
      auto itStart0 = itStart + (*itStart);
      auto itEnd0 = itStart + (*(itStart + 1));
      auto itStart1 = itStart + (*(itStart + 1));
      auto itEnd1 = itStart + (*(itStart + 2));
      auto itStart2 = itStart + (*(itStart + 2));
      auto itEnd2 = itStart + (*(itStart + 3));
      for(auto&& p: cooNew)
      {
         if(std::find(itStart0, itEnd0, p[0]) != itEnd0)
         {
            if(std::find(itStart1, itEnd1, p[1]) != itEnd1)
            {
               if(std::find(itStart2, itEnd2, p[2]) != itEnd2)
               {
                  cooFiltered.push_back(p);
                  count++;
               }
            }
         }
      }
      cooSizes.push_back(count);
   }
}

#if 1
void matchSendDispl(std::vector<std::vector<int>>& sendDispl, const std::vector<point_t>& locIdx, const std::vector<point_t>& remIdx, const std::vector<int>& remSizes, const std::vector<int>& remDispl)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::details::matchSendDispl");

   sendDispl.resize(remDispl.size());
   std::map<point_t, int> locOldIdx;
   for (std::size_t i = 0; i < locIdx.size(); ++i)
   {
      auto&& p = locIdx[i];
      locOldIdx[p] = i;
   }

   std::vector<point_t> locKeys;
   locKeys.reserve(locIdx.size());
   for(auto&& [p,v]: locOldIdx)
   {
      locKeys.emplace_back(p);
   }

   for (int r = 0; r < remDispl.size(); ++r)
   {
      // get new coo from other rank and check if it is here
      std::vector<point_t> remKeys;
      remKeys.reserve(remSizes.at(r));
      for (std::size_t i = 0; i < remSizes.at(r); ++i)
      {
         std::size_t ii = i + remDispl.at(r);
         auto&& p = remIdx[ii];
         remKeys.emplace_back(p);
      }
      std::sort(remKeys.begin(), remKeys.end());

      std::vector<point_t> sameKeys;
      sameKeys.reserve(std::min(remKeys.size(), locKeys.size()));
      if(remKeys.size() < locKeys.size())
      {
         std::set_intersection(remKeys.begin(), remKeys.end(), locKeys.begin(), locKeys.end(), std::back_inserter(sameKeys));
      }
      else
      {
         std::set_intersection(locKeys.begin(), locKeys.end(), remKeys.begin(), remKeys.end(), std::back_inserter(sameKeys));
      }

      sendDispl[r].reserve(sameKeys.size());
      for(auto&& p: sameKeys)
      {
         sendDispl[r].emplace_back(locOldIdx[p]);
      }
   }
}
#else
void matchSendDispl(std::vector<std::vector<int>>& sendDispl, const std::vector<point_t>& locIdx, const std::vector<point_t>& remIdx, const std::vector<int>& remSizes, const std::vector<int>& remDispl)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::details::matchSendDispl");

   sendDispl.resize(remDispl.size());
   std::map<point_t, int> locOldIdx;
   for (std::size_t i = 0; i < locIdx.size(); ++i)
   {
      auto&& p = locIdx[i];
      locOldIdx[p] = i;
   }

   for (int r = 0; r < remDispl.size(); ++r)
   {
      // get new coo from other rank and check if it is here
      std::map<point_t, int> remNewIdx;
      // setup remote map
      for (std::size_t i = 0; i < remSizes.at(r); ++i)
      {
         std::size_t ii = i + remDispl.at(r);
         auto&& p = remIdx[ii];
         remNewIdx[p] = i;
      }

      std::cerr << locOldIdx.size() << " vs " << remIdx.size() << std::endl;

      // loop over loc coo to find match
      for (auto itLCoo = locOldIdx.begin(); itLCoo != locOldIdx.end();)
      {
         auto lCoo = (*itLCoo).first;
         if (auto itRCoo = remNewIdx.find(lCoo); itRCoo != remNewIdx.end())
         {
            sendDispl[r].push_back((*itLCoo).second);
            itLCoo = locOldIdx.erase(itLCoo);
            remNewIdx.erase(itRCoo);
         }
         else
         {
            ++itLCoo;
         }
      }

      std::cerr << "OVERLAP: " << sendDispl[r].size() << std::endl;
   }
}
#endif
}

std::vector<std::vector<int>> getDispls(const std::vector<point_t>& absCooNew,
   const std::vector<point_t>& absCooOld, const MPI_Comm comm)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::getDispls");

   int rank, ranks;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &ranks);

   std::vector<int> locOldIdxSplit;
   details::getSplitIdx(locOldIdxSplit, absCooOld);

   std::vector<int> remNewIdxSplit;
   details::getSplitIdx(remNewIdxSplit, absCooNew);

   std::vector<int> remNewIdxSplitNeededAll;
   std::vector<int> remNewIdxSplitNeededSizes;
   for (int r = 0; r < ranks; ++r)
   {
      // get new coo from other rank and check if it is here
      std::vector<int> tmpSplit;
      const std::vector<int>* pRemSplit;

      pRemSplit = details::broadcastSplitIdx(r, rank, tmpSplit, remNewIdxSplit, comm);

      // Match split indexes
      details::matchSplitIdx(remNewIdxSplitNeededAll, remNewIdxSplitNeededSizes, locOldIdxSplit, *pRemSplit);
   }

   // Distributed split indexes
   std::vector<int> locNewIdxSplitNeededAll;
   std::vector<int> locNewIdxSplitNeededSizes;
   std::vector<int> locNewIdxSplitNeededDispl;
   details::distributeSplitIdx(remNewIdxSplitNeededAll, remNewIdxSplitNeededSizes, locNewIdxSplitNeededAll, locNewIdxSplitNeededSizes, locNewIdxSplitNeededDispl, comm);

   // Filter new indexes
   std::vector<point_t> absCooNewAll;
   std::vector<int> absCooNewSizes;
   details::filterIdx(absCooNewAll, absCooNewSizes, absCooNew, locNewIdxSplitNeededAll, locNewIdxSplitNeededDispl);

   // Distributed filtered split indexes
   std::vector<point_t> remAbsCooNewAll;
   std::vector<int> remAbsCooNewSizes;
   std::vector<int> remAbsCooNewDispl;
   details::distributeSplitIdx(absCooNewAll, absCooNewSizes, remAbsCooNewAll, remAbsCooNewSizes, remAbsCooNewDispl, comm);

   // Find matches for sendDispls
   std::vector<std::vector<int>> sendDispls;
   details::matchSendDispl(sendDispls, absCooOld, remAbsCooNewAll, remAbsCooNewSizes, remAbsCooNewDispl);

   return sendDispls;
}
#endif

std::vector<int> getReducedRanksSet(
   const std::vector<std::vector<int>>& sendDispls,
   const std::vector<std::vector<int>>& recvDispls, const MPI_Comm comm)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::getReducedRanksSet");

   std::set<int> commSet;
   std::set<int> sendSet;
   std::set<int> recvSet;
   // Save non-empty exchanges
   for (std::size_t i = 0; i < sendDispls.size(); ++i)
   {
      if (sendDispls[i].size() > 0)
      {
         commSet.insert(i);
         sendSet.insert(i);
      }
   }
   for (std::size_t i = 0; i < recvDispls.size(); ++i)
   {
      if (recvDispls[i].size() > 0)
      {
         commSet.insert(i);
         recvSet.insert(i);
      }
   }
   // Copy to vector
   std::vector<int> setLoc(commSet.size());
   std::size_t i = 0;
   for (auto it = commSet.begin(); it != commSet.end(); ++it)
   {
      setLoc[i++] = *it;
   }

   int rank, ranks;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &ranks);

   //
   // Check exchanges
   //

   // Recv remote set size
   std::vector<MPI_Request> req(sendSet.size() + recvSet.size());
   int count = 0;
   std::vector<int> remSetSize(sendSet.size(), 0);
   for (auto it = sendSet.begin(); it != sendSet.end(); ++it)
   {
      int sr = *it;
      MPI_Irecv(&remSetSize[count], 1, MPI_INT, sr, 0, comm, &req[count]);
      count++;
   }
   // Send size of remote set
   int size = setLoc.size();
   for (auto it = recvSet.begin(); it != recvSet.end(); ++it)
   {
      int rr = *it;
      MPI_Isend(&size, 1, MPI_INT, rr, 0, comm, &req[count]);
      count++;
   }
   // Wait for comm to be done
   std::vector<MPI_Status> stat(req.size());
   MPI_Waitall(req.size(), req.data(), stat.data());

   // Receive sets
   count = 0;
   std::vector<std::vector<int>> remSet(remSetSize.size());
   for (auto it = sendSet.begin(); it != sendSet.end(); ++it)
   {
      int sr = *it;
      remSet[count].resize(remSetSize[count]);
      MPI_Irecv(remSet[count].data(), remSet[count].size(), MPI_INT, sr, 1,
         comm, &req[count]);
      count++;
   }
   // Send local set
   for (auto it = recvSet.begin(); it != recvSet.end(); ++it)
   {
      int rr = *it;
      MPI_Isend(setLoc.data(), setLoc.size(), MPI_INT, rr, 1, comm,
         &req[count]);
      count++;
   }
   // Wait for comm to be done
   MPI_Waitall(req.size(), req.data(), stat.data());

   // Update local set
   for (auto& s: remSet)
   {
      for (auto r: s)
      {
         if (std::find(commSet.begin(), commSet.end(), r) == commSet.end())
         {
            // Remote set contained a rank missing form local set, add it
            commSet.insert(r);
         }
      }
   }

   // If local set was modified, update
   if (commSet.size() > setLoc.size())
   {
      setLoc.resize(commSet.size());
      std::size_t i = 0;
      for (auto it = commSet.begin(); it != commSet.end(); ++it)
      {
         setLoc[i++] = *it;
      }
   }

   return setLoc;
}

void redDisplsFromSet(std::vector<std::vector<int>>& sendDispls,
   std::vector<std::vector<int>>& recvDispls, const std::vector<int>& redSet)
{
   auto redSize = redSet.size();
   std::vector<std::vector<int>> sendDisplsRed(redSize);
   std::vector<std::vector<int>> recvDisplsRed(redSize);

   for (std::size_t r = 0; r < redSize; ++r)
   {
      sendDisplsRed[r] = std::move(sendDispls[redSet[r]]);
      recvDisplsRed[r] = std::move(recvDispls[redSet[r]]);
   }

   sendDispls = std::move(sendDisplsRed);
   recvDispls = std::move(recvDisplsRed);
}

MPI_Comm getSubComm(const std::vector<int>& redSet, const MPI_Comm comm)
{
   // Original group
   MPI_Group worldGroup;
   MPI_Comm_group(comm, &worldGroup);
   // Sub group
   MPI_Group subGroup;
   MPI_Group_incl(worldGroup, redSet.size(), redSet.data(), &subGroup);
   // Sub communicator
   MPI_Comm subComm;
   MPI_Comm_create(comm, subGroup, &subComm);
   return subComm;
}

} // namespace Mpi
} // namespace Transpose
} // namespace QuICC
