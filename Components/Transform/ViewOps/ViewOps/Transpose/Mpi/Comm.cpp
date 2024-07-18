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
#include "ViewOps/Transpose/Mpi/Comm.hpp"

namespace QuICC {
namespace Transpose {
namespace Mpi {


std::vector<std::vector<int>> getDispls(const std::vector<point_t>& absCooNew,
   const std::vector<point_t>& absCooOld, const MPI_Comm comm)
{
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
            // std::cout << "rank " << r << ": "
            //     << lCoo[0] << ',' << lCoo[1] << " found " << r << '\n';
            sendDispls[r].push_back((*itLCoo).second);
            itLCoo = locOldIdx.erase(itLCoo);
            remNewIdx.erase(itRCoo);
         }
         else
         {
            ++itLCoo;
            // std::cout << "not found\n";
         }
      }
   }
   return sendDispls;
}

std::vector<int> getReducedRanksSet(
   const std::vector<std::vector<int>>& sendDispls,
   const std::vector<std::vector<int>>& recvDispls, const MPI_Comm comm)
{
   std::set<int> redSet;
   // Save non-empty exchanges
   for (std::size_t i = 0; i < sendDispls.size(); ++i)
   {
      if (sendDispls[i].size() > 0)
      {
         redSet.insert(i);
      }
   }
   for (std::size_t i = 0; i < recvDispls.size(); ++i)
   {
      if (recvDispls[i].size() > 0)
      {
         redSet.insert(i);
      }
   }
   // Copy to vector
   std::vector<int> setLoc(redSet.size());
   std::size_t i = 0;
   for (auto it = redSet.begin(); it != redSet.end(); ++it)
   {
      setLoc[i++] = *it;
   }
   // Get max set size and location
   int rank;
   MPI_Comm_rank(comm, &rank);
   std::array<int, 2> sizeLoc = {static_cast<int>(setLoc.size()), rank};
   std::array<int, 2> sizeLocGlo;
   MPI_Allreduce(sizeLoc.data(), sizeLocGlo.data(), 1, MPI_2INT, MPI_MAXLOC,
      comm);
   // Get global set
   std::vector<int> setGlo(sizeLocGlo[0]);
   int rankMax = sizeLocGlo[1];
   if (rank == rankMax)
   {
      MPI_Bcast(setLoc.data(), setLoc.size(), MPI_INT, rankMax, comm);
      return setLoc;
   }
   else
   {
      MPI_Bcast(setGlo.data(), setGlo.size(), MPI_INT, rankMax, comm);
      // This rank is not in the set
      if (setLoc.size() == 0)
      {
         return setLoc;
      }
      // Check set intersection
      for (auto l : setLoc)
      {
         assert(std::find(setGlo.begin(), setGlo.end(), l) != setGlo.end());
      }
      return setGlo;
   }
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

std::vector<int> getCount(const std::vector<std::vector<int>>& displs)
{
   std::vector<int> count(displs.size());
   for (std::size_t i = 0; i < displs.size(); ++i)
   {
      if (displs[i].size() > 0)
      {
         count[i] = 1;
      }
      else
      {
         count[i] = 0;
      }
   }
   return count;
}

} // namespace Mpi
} // namespace Transpose
} // namespace QuICC
