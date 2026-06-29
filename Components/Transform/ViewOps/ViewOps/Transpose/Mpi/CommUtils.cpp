/**
 * @file Comm.cpp
 * @brief Methods for mpi enabled transform
 */

// External includes
//
#include <algorithm>
#include <numeric>
#include <set>

// Project includes
//
#include "ViewOps/Transpose/Mpi/CommUtils.hpp"
#include "Utils/Utils.hpp"
#include "Utils/Mpi/Utils.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Transpose {
namespace Mpi {

std::vector<std::vector<int>> getDispls(const std::vector<point_t>& absCooNew,
   const std::vector<point_t>& absCooOld, const MPI_Comm comm)
{
   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::getDispls");

   int rank, ranks;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &ranks);

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

      pRemSplit = Utils::Mpi::broadcastSplitIdx(r, rank, tmpSplit, remNewIdxSplit, comm);

      // Match split indexes
      Utils::matchSplitIdx(remNewIdxSplitNeededAll, remNewIdxSplitNeededSizes, locOldIdxSplit, *pRemSplit);
   }

   // Distributed split indexes
   std::vector<int> locNewIdxSplitNeededAll;
   std::vector<int> locNewIdxSplitNeededSizes;
   std::vector<int> locNewIdxSplitNeededDispl;
   Utils::Mpi::distributeSplitIdx(remNewIdxSplitNeededAll, remNewIdxSplitNeededSizes, locNewIdxSplitNeededAll, locNewIdxSplitNeededSizes, locNewIdxSplitNeededDispl, comm);

   // Filter new indexes
   std::vector<point_t> absCooNewAll;
   std::vector<int> absCooNewSizes;
   Utils::filterIdx(absCooNewAll, absCooNewSizes, absCooNew, locNewIdxSplitNeededAll, locNewIdxSplitNeededDispl);

   // Distributed filtered split indexes
   std::vector<point_t> remAbsCooNewAll;
   std::vector<int> remAbsCooNewSizes;
   std::vector<int> remAbsCooNewDispl;
   Utils::Mpi::distributeSplitIdx(absCooNewAll, absCooNewSizes, remAbsCooNewAll, remAbsCooNewSizes, remAbsCooNewDispl, comm);

   // Find matches for sendDispls
   std::vector<std::vector<int>> sendDispls;
   Utils::matchSendDispl(sendDispls, absCooOld, remAbsCooNewAll, remAbsCooNewSizes, remAbsCooNewDispl);

   return sendDispls;
}

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
   // Check exchanges (iterate to catch multiple hop connections)
   //

   int commSet_changed = 1;
   while(commSet_changed)
   {
      // Traverse graph (catches nodes which only send)
      commSet_changed = collectRemoteSet(setLoc, commSet, sendSet, recvSet, comm);

      // Traverse graph in reverse direction (catches nodes which only receive)
      commSet_changed += collectRemoteSet(setLoc, commSet, recvSet, sendSet, comm);

      // Global changes
      MPI_Allreduce(MPI_IN_PLACE, &commSet_changed, 1, MPI_INT, MPI_SUM, comm);
   }

   return setLoc;
}

bool collectRemoteSet(std::vector<int>& setLoc, std::set<int>& commSet, const std::set<int>& sendSet, const std::set<int>& recvSet, const MPI_Comm comm)
{
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
   bool changedSet = false;
   if (commSet.size() > setLoc.size())
   {
      changedSet = true;

      setLoc.resize(commSet.size());
      std::size_t i = 0;
      for (auto it = commSet.begin(); it != commSet.end(); ++it)
      {
         setLoc[i++] = *it;
      }
   }

   return changedSet;
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
