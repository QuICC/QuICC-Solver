#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>
#include <chrono>
#include <climits>
#include <regex>
#include <thread>

#include "ViewOps/Transpose/MpiUtils.hpp"
extern "C" {
#include <unistd.h>
}


using namespace QuICC::Transpose;


int main(int argc, char** argv)
{
   MPI_Init(NULL, NULL);

   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   // if(rank == 0)
   // {
   //     using namespace std::chrono_literals;
   //     volatile bool wait = true;
   //     char hostname[HOST_NAME_MAX];
   //     gethostname(hostname, HOST_NAME_MAX);
   //     std::cerr << "PID " << getpid() << " on " << hostname
   //         << " ready to attach" << std::endl;
   //     while (wait == true)
   //     std::this_thread::sleep_for(1s);
   // }
   // MPI_Barrier(MPI_COMM_WORLD);

   Catch::Session session; // There must be exactly one instance

   auto returnCode = session.run();

   MPI_Finalize();

   return returnCode;
}


TEST_CASE("Send/recv displacements", "[SendRecvDispls]")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   assert(ranks >= 2);

   // this should be store in a class
   std::vector<int> sendBuf;
   std::vector<int> recvBuf;
   std::vector<int> sendCounts(ranks, 1);
   std::vector<int> recvCounts(ranks, 1);
   std::vector<int> sDispls(ranks, 0);
   std::vector<int> rDispls(ranks, 0);
   std::vector<std::vector<int>> sendDisplsRef(ranks);
   std::vector<std::vector<int>> recvDisplsRef(ranks);

   std::vector<MPI_Datatype> sendType(ranks);
   std::vector<MPI_Datatype> recvType(ranks);

   // needed to build sendRankMap recvRankMap
   std::vector<point_t> absCooOld;
   std::vector<point_t> absCooNew;


   //
   // transpose distributed trapezoidal matrix

   // r0   0  1  *  *    r1   *  *  2  3
   //      4  5  *  *         *  *  6  7
   //      8  9  *            *  * 10

   //             vvv to vvv

   // r0   0  4  *        r1  *  *  8
   //      1  5  *            *  *  9
   //      2  6  *            *  * 10
   //      3  7

   if (rank == 0)
   {
      absCooOld = {{0, 0, 0}, {0, 1, 0}, {0, 2, 0}, {1, 0, 0}, {1, 1, 0},
         {1, 2, 0}};
      absCooNew = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {0, 1, 0},
         {1, 1, 0}, {2, 1, 0}, {3, 1, 0}};
      sendBuf = {0, 4, 8, 1, 5, 9}; // col maj
      recvBuf = std::vector<int>(8, -1);
      sendDisplsRef[0] = {0, 1, 3, 4}; //
      recvDisplsRef[0] = {0, 4, 1, 5}; // 0 4 1 5
      sendDisplsRef[1] = {2, 5};       //
      recvDisplsRef[1] = {2, 6, 3, 7}; // 2 6 3 7
   }
   if (rank == 1)
   {
      absCooOld = {{2, 0, 0}, {2, 1, 0}, {2, 2, 0}, {3, 0, 0}, {3, 1, 0}};
      absCooNew = {{0, 2, 0}, {1, 2, 0}, {2, 2, 0}};
      sendBuf = {2, 6, 10, 3, 7}; // col maj
      recvBuf = std::vector<int>(3, -1);
      sendDisplsRef[0] = {0, 1, 3, 4}; //
      recvDisplsRef[0] = {0, 1};       // 8 9
      sendDisplsRef[1] = {2};          //
      recvDisplsRef[1] = {2};          // 10
   }

   auto sendDispls = getDispls(absCooNew, absCooOld, MPI_COMM_WORLD);
   auto recvDispls = getDispls(absCooOld, absCooNew, MPI_COMM_WORLD);

   // Check
   for (int r = 0; r < ranks; ++r)
   {
      CHECK(sendDispls[r].size() == sendDisplsRef[r].size());
      for (std::size_t i = 0; i < sendDispls[r].size(); ++i)
      {
         CHECK(sendDispls[r][i] == sendDisplsRef[r][i]);
      }
   }
   for (int r = 0; r < ranks; ++r)
   {
      CHECK(recvDispls[r].size() == recvDisplsRef[r].size());
      for (std::size_t i = 0; i < recvDispls[r].size(); ++i)
      {
         CHECK(recvDispls[r][i] == recvDisplsRef[r][i]);
      }
   }
}

TEST_CASE("Reduce displacements", "[ReduceDispls]")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   std::vector<std::vector<int>> sendDispls(ranks);
   std::vector<std::vector<int>> recvDispls(ranks);
   std::vector<std::vector<int>> sendDisplsRef(ranks);
   std::vector<std::vector<int>> recvDisplsRef(ranks);

   if (rank == 0)
   {
      sendDispls[0] = {0, 1, 3, 4};
      recvDispls[0] = {0, 4, 1, 5};
      sendDispls[1] = {2, 5};
      recvDispls[1] = {2, 6, 3, 7};
      sendDisplsRef[0] = {0, 1, 3, 4};
      recvDisplsRef[0] = {0, 4, 1, 5};
      sendDisplsRef[1] = {2, 5};
      recvDisplsRef[1] = {2, 6, 3, 7};
   }
   if (rank == 1)
   {
      sendDispls[0] = {0, 1, 3, 4};
      recvDispls[0] = {0, 1};
      sendDispls[1] = {2};
      recvDispls[1] = {2};
      sendDisplsRef[0] = {0, 1, 3, 4};
      recvDisplsRef[0] = {0, 1};
      sendDisplsRef[1] = {2};
      recvDisplsRef[1] = {2};
   }

   auto redSet = getReducedRanksSet(sendDispls, recvDispls);

   if (rank == 0 || rank == 1)
   {
      CHECK(redSet.size() == 2);
   }
   else
   {
      CHECK(redSet.size() == 0);
   }

   redDisplsFromSet(sendDispls, recvDispls, redSet);

   if (rank == 0 || rank == 1)
   {
      CHECK(sendDispls.size() == 2);
      CHECK(recvDispls.size() == 2);
   }
   else
   {
      CHECK(sendDispls.size() == 0);
      CHECK(recvDispls.size() == 0);
   }

   for (std::size_t r = 0; r < redSet.size(); ++r)
   {
      for (std::size_t i = 0; i < sendDispls[r].size(); ++i)
      {
         CHECK(sendDispls[r][i] == sendDisplsRef[redSet[r]][i]);
      }
      for (std::size_t i = 0; i < recvDispls[r].size(); ++i)
      {
         CHECK(recvDispls[r][i] == recvDisplsRef[redSet[r]][i]);
      }
   }
}


TEST_CASE("SubComm and Types", "[SubComm]")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   assert(ranks >= 3);

   // this should be store in a class
   std::vector<int> sendBuf;
   std::vector<int> recvBuf;
   std::vector<int> recvBufRef;
   std::vector<std::vector<int>> sendDisplsRef(ranks);
   std::vector<std::vector<int>> recvDisplsRef(ranks);

   // needed to build sendRankMap recvRankMap
   std::vector<point_t> absCooOld;
   std::vector<point_t> absCooNew;

   //
   // transpose distributed trapezoidal matrix
   // with one sided comm

   // r0   0  1  *  *  *   r1   *  *  2  3   r2 * * * * 4
   //      5  6  *  *           *  *  7  8      * * * *
   //      9 10  *              *  * 11         * * *

   //             vvv to vvv

   // r0   0  5  *        r1  *  *  9        r2  * * *
   //      1  6  *            *  * 10            * * *
   //      2  7  *            *  * 11            * * *
   //      3  8               *  *               * *
   //      4                  *                  *

   if (rank == 0)
   {
      absCooOld = {{0, 0, 0}, {0, 1, 0}, {0, 2, 0}, {1, 0, 0}, {1, 1, 0},
         {1, 2, 0}};
      absCooNew = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {4, 0, 0},
         {0, 1, 0}, {1, 1, 0}, {2, 1, 0}, {3, 1, 0}};
      sendBuf = {0, 5, 9, 1, 6, 10}; // col maj
      recvBuf = std::vector<int>(9, -1);
      recvBufRef = {0, 1, 2, 3, 4, 5, 6, 7, 8};
   }
   if (rank == 1)
   {
      absCooOld = {{2, 0, 0}, {2, 1, 0}, {2, 2, 0}, {3, 0, 0}, {3, 1, 0}};
      absCooNew = {{0, 2, 0}, {1, 2, 0}, {2, 2, 0}};
      sendBuf = {2, 7, 11, 3, 8}; // col maj
      recvBuf = std::vector<int>(3, -1);
      recvBufRef = {9, 10, 11};
   }
   if (rank == 2)
   {
      absCooOld = {{4, 0, 0}};
      sendBuf = {4}; // col maj
   }

   auto sendDispls = getDispls(absCooNew, absCooOld);
   auto recvDispls = getDispls(absCooOld, absCooNew);
   auto redSet = getReducedRanksSet(sendDispls, recvDispls);

   // Check
   if (rank <= 2)
   {
      CHECK(redSet.size() == 3);
   }
   else
   {
      CHECK(redSet.size() == 0);
   }

   redDisplsFromSet(sendDispls, recvDispls, redSet);
   auto subComm = getSubComm(redSet);

   if (rank <= 2)
   {
      CHECK(subComm != MPI_COMM_NULL);
   }
   else
   {
      CHECK(subComm == MPI_COMM_NULL);
   }

   // Communicate with partecipants only
   if (subComm != MPI_COMM_NULL)
   {

      // Setup comm
      int subRanks;
      MPI_Comm_size(subComm, &subRanks);

      std::vector<MPI_Datatype> sendType(subRanks);
      std::vector<MPI_Datatype> recvType(subRanks);

      // Build types
      for (int r = 0; r < subRanks; ++r)
      {
         MPI_Type_create_indexed_block(sendDispls[r].size(), 1,
            sendDispls[r].data(), MPI_INT, &sendType[r]);
         MPI_Type_commit(&sendType[r]);
         MPI_Type_create_indexed_block(recvDispls[r].size(), 1,
            recvDispls[r].data(), MPI_INT, &recvType[r]);
         MPI_Type_commit(&recvType[r]);
      }

      auto sendCounts = getCount(sendDispls);
      auto recvCounts = getCount(recvDispls);

      std::vector<int> sDispls(subRanks, 0);
      std::vector<int> rDispls(subRanks, 0);

      // Comm
      MPI_Alltoallw(sendBuf.data(), sendCounts.data(), sDispls.data(),
         sendType.data(), recvBuf.data(), recvCounts.data(), rDispls.data(),
         recvType.data(), subComm);

      // Check
      for (std::size_t i = 0; i < recvBufRef.size(); ++i)
      {
         CHECK(recvBuf[i] == recvBufRef[i]);
      }

      // Release types
      for (int r = 0; r < subRanks; ++r)
      {
         MPI_Type_free(&sendType[r]);
         MPI_Type_free(&recvType[r]);
      }
   }
}
