/**
 * @file Utils.cpp
 * @brief General MPI utils
 */

// External includes
//
#include <algorithm>
#include <numeric>
#include <set>

// Project includes
//
#include "Utils/Mpi/Utils.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Utils {
namespace Mpi {

const std::vector<int>* broadcastSplitIdx(const int r, const int rank, std::vector<int>& remSplit, const std::vector<int>& locSplit, const MPI_Comm comm)
{
   Profiler::RegionFixture<4> fix("Utils::broadcastSplitIdx");

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
   Profiler::RegionFixture<5> fix("Utils::distributeSplitSizes");

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
   Profiler::RegionFixture<4> fix("Utils::distributeSplitIdx");

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
   Profiler::RegionFixture<4> fix("Utils::distributeSplitIdx_point");

   const auto dimSize = std::tuple_size<point_t>{};

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

} // namespace Mpi
} // namespace Utils
} // namespace QuICC
